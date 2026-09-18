function actualPath = exportQSSResults(result,outputPath,opts)
%EXPORTQSSRESULTS Atomically write a validated QSS telemetry v1 HDF5 file.
%   Existing output names are made unique by default. Set opts.overwrite=true
%   to replace the requested path explicitly.

if nargin < 3 || isempty(opts), opts = struct(); end
if ~isstruct(opts) || numel(opts) ~= 1
    error('exportQSSResults:badOptions','opts must be a scalar struct.');
end
if ~isfield(opts,'overwrite') || isempty(opts.overwrite), opts.overwrite = false; end
if ~isfield(opts,'unique_name') || isempty(opts.unique_name), opts.unique_name = true; end
opts.overwrite = logical(opts.overwrite); opts.unique_name = logical(opts.unique_name);
[ok,report] = validateQSSResults(result);
if ~ok
    error('exportQSSResults:invalidResult','Result failed validation: %s',strjoin(report.errors,'; '));
end
target = resolveTarget(outputPath,opts);
folder = fileparts(target);
if isempty(folder), folder = pwd; end
if ~isfolder(folder), mkdir(folder); end
[~,base,ext] = fileparts(target);
if isempty(ext), ext = '.h5'; target = [target ext]; end
tmp = fullfile(folder,['.' base '.qss-' temporaryId() '.tmp' ext]);
cleanup = onCleanup(@() deleteIfPresent(tmp)); %#ok<NASGU>
try
    writeHdf5(result,tmp);
    validateWrittenFile(tmp,result);
catch errorIn
    error('exportQSSResults:writeFailed','Temporary file validation failed: %s',errorIn.message);
end
[moved,message] = movefile(tmp,target,'f');
if ~moved
    error('exportQSSResults:replaceFailed','Could not replace output %s: %s',target,message);
end
actualPath = char(target);
end

function target = resolveTarget(outputPath,opts)
if ~(ischar(outputPath) || (isstring(outputPath) && isscalar(outputPath))) || strlength(string(outputPath))==0
    error('exportQSSResults:badPath','outputPath must be a non-empty path.');
end
target = char(string(outputPath));
[folder,base,ext] = fileparts(target);
if isempty(ext), ext = '.h5'; target = fullfile(folder,[base ext]); end
if isfile(target) && ~opts.overwrite
    if ~opts.unique_name
        error('exportQSSResults:exists','Output exists; set overwrite=true or unique_name=true.');
    end
    for i = 1:100000
        candidate = fullfile(folder,sprintf('%s_%03d%s',base,i,ext));
        if ~isfile(candidate), target = candidate; return; end
    end
    error('exportQSSResults:noUniqueName','Could not choose a unique output name.');
end
end

function writeHdf5(result,path)
writeJson(path,'/manifest_json',result.manifest);
for i = 1:numel(result.tracks)
    t = result.tracks(i); root = ['/tracks/'h5Name(t.id)];
    writeJson(path,[root '/metadata_json'],t.metadata);
    writeVector(path,[root '/distance_m'],t.distance_m,'double');
    writeVector(path,[root '/curvature_per_m'],t.curvature_per_m,'double');
    if isfield(t,'x_m') && ~isempty(t.x_m), writeVector(path,[root '/x_m'],t.x_m,'double'); end
    if isfield(t,'y_m') && ~isempty(t.y_m), writeVector(path,[root '/y_m'],t.y_m,'double'); end
end
for i = 1:numel(result.cases)
    c = result.cases(i); croot = ['/cases/'h5Name(c.id)];
    writeJson(path,[croot '/metadata_json'],c.metadata);
    writeJson(path,[croot '/setup_json'],c.setup);
    if isfield(c,'envelope') && ~isempty(c.envelope)
        for j = 1:numel(c.envelope)
            writeEnvelope(path,[croot '/envelope/' h5Name(c.envelope(j).id)],c.envelope(j));
        end
    end
    for j = 1:numel(c.laps)
        writeLap(path,[croot '/laps/' h5Name(c.laps(j).id)],c.laps(j));
    end
end
end

function writeLap(path,root,lap)
writeJson(path,[root '/metadata_json'],lap.metadata);
for i = 1:numel(lap.axes)
    a = lap.axes(i); ar = [root '/axes/' h5Name(a.id)];
    writeVector(path,[ar '/time_s'],a.time_s,'double');
    writeVector(path,[ar '/distance_m'],a.distance_m,'double');
end
for i = 1:numel(lap.channels)
    writeChannel(path,[root '/channels/' h5Name(lap.channels(i).id)],lap.channels(i));
end
if isfield(lap,'diagnostics') && ~isempty(lap.diagnostics)
    writeJson(path,[root '/diagnostics_json'],lap.diagnostics);
end
if isfield(lap,'reconstruction') && ~isempty(lap.reconstruction)
    writeJson(path,[root '/reconstruction_json'],lap.reconstruction);
end
end

function writeEnvelope(path,root,ch)
writeVector(path,[root '/values'],ch.values,'double');
writeVector(path,[root '/valid'],uint8(logical(ch.valid)),'uint8');
writeJson(path,[root '/metadata_json'],channelMetadata(ch));
end

function writeChannel(path,root,ch)
writeVector(path,[root '/values'],ch.values,'double');
writeVector(path,[root '/valid'],uint8(logical(ch.valid)),'uint8');
writeJson(path,[root '/metadata_json'],channelMetadata(ch));
end

function metadata = channelMetadata(ch)
metadata = struct('id',char(string(ch.id)),'label',char(string(ch.label)), ...
    'unit',char(string(ch.unit)),'axis_id',char(string(ch.axis_id)), ...
    'origin',char(string(ch.origin)),'interpolation',char(string(ch.interpolation)), ...
    'description',char(string(ch.description)),'coordinate_frame',char(string(ch.coordinate_frame)), ...
    'sign_convention',char(string(ch.sign_convention)));
end

function writeJson(path,dataset,value)
bytes = jsonBytes(value);
h5create(path,dataset,numel(bytes),'Datatype','uint8');
h5write(path,dataset,bytes);
end

function writeVector(path,dataset,value,datatype)
value = value(:);
% Numeric telemetry arrays are chunked and compressed as required by the
% interchange contract.  Keep chunks bounded so small fixture files and
% long reconstructed laps use the same portable layout.
chunkSize = max(1,min(numel(value),1024));
h5create(path,dataset,numel(value),'Datatype',datatype,'ChunkSize',chunkSize,'Deflate',6);
h5write(path,dataset,value);
end

function bytes = jsonBytes(value)
value = safeJson(value);
textValue = jsonencode(value);
bytes = uint8(unicode2native(textValue,'UTF-8'));
bytes = bytes(:);
end

function value = safeJson(input)
if isstruct(input)
    value = repmat(struct(),size(input)); fields = fieldnames(input);
    for n = 1:numel(input)
        for i = 1:numel(fields)
            try, value(n).(fields{i}) = safeJson(input(n).(fields{i}));
            catch, value(n).(fields{i}) = class(input(n).(fields{i}));
            end
        end
    end
elseif iscell(input)
    value = cell(size(input)); for i=1:numel(input), value{i}=safeJson(input{i}); end
elseif isobject(input)
    value = class(input);
elseif isa(input,'function_handle')
    value = func2str(input);
elseif isstring(input)
    value = cellstr(input);
elseif isnumeric(input)
    value = double(input); value(~isfinite(value)) = NaN;
elseif ischar(input) || islogical(input)
    value = input;
else
    value = char(string(input));
end
end

function validateWrittenFile(path,result)
try
    manifest = readJson(path,'/manifest_json');
catch err
    error('exportQSSResults:manifestRead','Could not read UTF-8 manifest: %s',err.message);
end
if ~isstruct(manifest) || ~isfield(manifest,'format') || ~strcmp(char(string(manifest.format)),'qss-telemetry')
    error('exportQSSResults:badManifest','Written manifest is not qss-telemetry.');
end
if ~isfield(manifest,'schema_version') || ~strcmp(char(string(manifest.schema_version)),'1.0')
    error('exportQSSResults:badManifestVersion','Written schema version is not 1.0.');
end
for i = 1:numel(result.tracks)
    root = ['/tracks/'h5Name(result.tracks(i).id)];
    mustExist(path,[root '/metadata_json']); mustExist(path,[root '/distance_m']); mustExist(path,[root '/curvature_per_m']);
end
for i = 1:numel(result.cases)
    c = result.cases(i); root = ['/cases/'h5Name(c.id)];
    mustExist(path,[root '/metadata_json']); mustExist(path,[root '/setup_json']);
    for j = 1:numel(c.laps)
        lap = c.laps(j); lr = [root '/laps/'h5Name(lap.id)];
        mustExist(path,[lr '/metadata_json']);
        for a = 1:numel(lap.axes)
            ar = [lr '/axes/' h5Name(lap.axes(a).id)];
            mustExist(path,[ar '/time_s']); mustExist(path,[ar '/distance_m']);
        end
        for k = 1:numel(lap.channels)
            cr = [lr '/channels/' h5Name(lap.channels(k).id)];
            mustExist(path,[cr '/values']); mustExist(path,[cr '/valid']); mustExist(path,[cr '/metadata_json']);
            values = h5read(path,[cr '/values']); valid = h5read(path,[cr '/valid']);
            if numel(values) ~= numel(lap.channels(k).values) || numel(valid) ~= numel(values)
                error('exportQSSResults:shapeMismatch','Written channel shape mismatch at %s.',cr);
            end
        end
    end
end
end

function value = readJson(path,dataset)
bytes = uint8(h5read(path,dataset)); bytes = bytes(:).';
value = jsondecode(native2unicode(bytes,'UTF-8'));
end

function mustExist(path,dataset)
try
    h5info(path,dataset);
catch err
    error('exportQSSResults:missingDataset','Missing dataset %s: %s',dataset,err.message);
end
end

function name = h5Name(value)
name = char(string(value));
if isempty(name) || contains(name,'/') || any(name==0)
    error('exportQSSResults:unsafeId','IDs used as HDF5 paths may not be empty or contain /.');
end
end

function value = temporaryId()
try, value = char(java.util.UUID.randomUUID()); catch, value = sprintf('%d',randi(2^31-1)); end
value = strrep(value,'-','');
end

function deleteIfPresent(path)
if isfile(path), delete(path); end
end
