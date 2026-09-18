function [isValid,report] = validateQSSResults(result)
%VALIDATEQSSRESULTS Validate the in-memory QSS telemetry v1 interchange.
%   [OK,REPORT] = validateQSSResults(RESULT) checks the MATLAB struct used by
%   buildTelemetryResult/exportQSSResults. With no output arguments an invalid
%   result raises an error, which is convenient at export boundaries while
%   still allowing tests and callers to inspect all failures.

errors = {};
warnings = {};

if ~isstruct(result) || numel(result) ~= 1
    errors{end+1} = 'result must be a scalar struct'; %#ok<AGROW>
else
    errors = checkManifest(result,errors);
    if isfield(result,'tracks')
        [errors,warnings] = checkTracks(result.tracks,errors,warnings);
    else
        errors{end+1} = 'result.tracks is required'; %#ok<AGROW>
    end
    if isfield(result,'cases')
        [errors,warnings] = checkCases(result.cases,result,errors,warnings);
    else
        errors{end+1} = 'result.cases is required'; %#ok<AGROW>
    end
end

isValid = isempty(errors);
report = struct('valid',isValid,'errors',{errors},'warnings',{warnings});
if nargout == 0 && ~isValid
    error('validateQSSResults:invalidResult','Invalid QSS telemetry result: %s', ...
        strjoin(errors,'; '));
end
end

function errors = checkManifest(result,errors)
if ~isfield(result,'manifest') || ~isstruct(result.manifest)
    errors{end+1} = 'manifest is required'; %#ok<AGROW>
    return
end
m = result.manifest;
required = {'format','schema_version','run_uuid','creation_utc','source_type', ...
    'producer','matlab_version','git_revision','git_dirty','tracks','cases'};
for i = 1:numel(required)
    if ~isfield(m,required{i})
        errors{end+1} = sprintf('manifest.%s is required',required{i}); %#ok<AGROW>
    end
end
if isfield(m,'format') && ~strcmp(char(string(m.format)),'qss-telemetry')
    errors{end+1} = 'manifest.format must be qss-telemetry'; %#ok<AGROW>
end
if isfield(m,'schema_version') && ~strcmp(char(string(m.schema_version)),'1.0')
    errors{end+1} = 'manifest.schema_version must be 1.0'; %#ok<AGROW>
end
end

function [errors,warnings] = checkTracks(tracks,errors,warnings)
if ~isstruct(tracks) || isempty(tracks)
    errors{end+1} = 'tracks must be a non-empty struct array'; %#ok<AGROW>
    return
end
ids = strings(0,1);
for i = 1:numel(tracks)
    t = tracks(i);
    prefix = sprintf('tracks(%d)',i);
    [ok,id] = textField(t,'id');
    if ~ok || strlength(id) == 0 || contains(id,'/')
        errors{end+1} = [prefix '.id must be a non-empty path-safe string']; %#ok<AGROW>
    else
        if any(ids == id)
            errors{end+1} = sprintf('duplicate track id %s',id); %#ok<AGROW>
        end
        ids(end+1,1) = id; %#ok<AGROW>
    end
    [errors,~,~] = checkFiniteMonotonicVector(t,prefix,'distance_m',errors,true);
    [errors,~,~] = checkFiniteMonotonicVector(t,prefix,'curvature_per_m',errors,false);
    if isfield(t,'distance_m') && isfield(t,'curvature_per_m') && ...
            isnumeric(t.distance_m) && isnumeric(t.curvature_per_m) && ...
            numel(t.distance_m) ~= numel(t.curvature_per_m)
        errors{end+1} = [prefix '.distance_m and curvature_per_m must have equal length']; %#ok<AGROW>
    end
    hasX = isfield(t,'x_m') && ~isempty(t.x_m);
    hasY = isfield(t,'y_m') && ~isempty(t.y_m);
    if xor(hasX,hasY)
        errors{end+1} = [prefix '.x_m and y_m must either both exist or both be absent']; %#ok<AGROW>
    elseif hasX
        if ~isnumeric(t.x_m) || ~isnumeric(t.y_m) || ...
                numel(t.x_m) ~= numel(t.distance_m) || numel(t.y_m) ~= numel(t.distance_m)
            errors{end+1} = [prefix '.x_m and y_m must match distance_m']; %#ok<AGROW>
        end
    end
end
end

function [errors,warnings] = checkCases(cases,result,errors,warnings)
if ~isstruct(cases) || isempty(cases)
    errors{end+1} = 'cases must be a non-empty struct array'; %#ok<AGROW>
    return
end
trackIds = strings(0,1);
if isfield(result,'tracks') && isstruct(result.tracks)
    for i = 1:numel(result.tracks)
        if isfield(result.tracks(i),'id')
            trackIds(end+1,1) = string(result.tracks(i).id); %#ok<AGROW>
        end
    end
end
caseIds = strings(0,1);
for i = 1:numel(cases)
    c = cases(i);
    prefix = sprintf('cases(%d)',i);
    [ok,id] = textField(c,'id');
    if ~ok || strlength(id) == 0 || contains(id,'/')
        errors{end+1} = [prefix '.id must be a non-empty path-safe string']; %#ok<AGROW>
    else
        if any(caseIds == id)
            errors{end+1} = sprintf('duplicate case id %s',id); %#ok<AGROW>
        end
        caseIds(end+1,1) = id; %#ok<AGROW>
    end
    if ~isfield(c,'laps') || ~isstruct(c.laps) || isempty(c.laps)
        errors{end+1} = [prefix '.laps must be a non-empty struct array']; %#ok<AGROW>
        continue
    end
    lapIds = strings(0,1);
    for j = 1:numel(c.laps)
        lap = c.laps(j);
        lprefix = sprintf('%s.laps(%d)',prefix,j);
        [ok,lapId] = textField(lap,'id');
        if ~ok || strlength(lapId) == 0 || contains(lapId,'/')
            errors{end+1} = [lprefix '.id must be a non-empty path-safe string']; %#ok<AGROW>
        elseif any(lapIds == lapId)
            errors{end+1} = sprintf('duplicate lap id %s in case %s',lapId,id); %#ok<AGROW>
        else
            lapIds(end+1,1) = lapId; %#ok<AGROW>
        end
        if isfield(lap,'track_id') && ~isempty(lap.track_id) && ...
                ~any(trackIds == string(lap.track_id))
            errors{end+1} = sprintf('%s.track_id refers to missing track',lprefix); %#ok<AGROW>
        end
        [errors,warnings] = checkLap(lap,lprefix,errors,warnings);
    end
end
end

function [errors,warnings] = checkLap(lap,prefix,errors,warnings)
if ~isfield(lap,'axes') || ~isstruct(lap.axes) || isempty(lap.axes)
    errors{end+1} = [prefix '.axes must be a non-empty struct array']; %#ok<AGROW>
    return
end
axisIds = strings(0,1);
for i = 1:numel(lap.axes)
    a = lap.axes(i);
    aprefix = sprintf('%s.axes(%d)',prefix,i);
    [ok,id] = textField(a,'id');
    if ~ok || strlength(id) == 0 || contains(id,'/')
        errors{end+1} = [aprefix '.id must be a non-empty path-safe string']; %#ok<AGROW>
    elseif any(axisIds == id)
        errors{end+1} = sprintf('duplicate axis id %s',id); %#ok<AGROW>
    else
        axisIds(end+1,1) = id; %#ok<AGROW>
    end
    [errors,timeSize,~] = checkFiniteMonotonicVector(a,aprefix,'time_s',errors,true);
    [errors,distSize,~] = checkFiniteMonotonicVector(a,aprefix,'distance_m',errors,true);
    if ~isempty(timeSize) && ~isempty(distSize) && timeSize ~= distSize
        errors{end+1} = [aprefix '.time_s and distance_m must have equal length']; %#ok<AGROW>
    end
end
if ~isfield(lap,'channels') || isempty(lap.channels)
    warnings{end+1} = [prefix '.channels is empty']; %#ok<AGROW>
    return
end
if ~isstruct(lap.channels)
    errors{end+1} = [prefix '.channels must be a struct array']; %#ok<AGROW>
    return
end
channelIds = strings(0,1);
for i = 1:numel(lap.channels)
    ch = lap.channels(i);
    cprefix = sprintf('%s.channels(%d)',prefix,i);
    [ok,id] = textField(ch,'id');
    if ~ok && isfield(ch,'metadata') && isstruct(ch.metadata)
        [ok,id] = textField(ch.metadata,'id');
    end
    if ~ok || strlength(id) == 0 || contains(id,'/')
        errors{end+1} = [cprefix '.id must be a non-empty path-safe string']; %#ok<AGROW>
    elseif any(channelIds == id)
        errors{end+1} = sprintf('duplicate channel id %s',id); %#ok<AGROW>
    else
        channelIds(end+1,1) = id; %#ok<AGROW>
    end
    metadata = ch;
    if isfield(ch,'metadata') && isstruct(ch.metadata), metadata = ch.metadata; end
    required = {'id','label','unit','axis_id','origin','interpolation', ...
        'description','coordinate_frame','sign_convention'};
    for k = 1:numel(required)
        if ~isfield(metadata,required{k})
            errors{end+1} = sprintf('%s metadata.%s is required',cprefix,required{k}); %#ok<AGROW>
        end
    end
    if isfield(metadata,'origin') && ~any(strcmp(char(string(metadata.origin)), ...
            {'simulation_output','derived','qss_reconstructed','measured'}))
        errors{end+1} = sprintf('%s has invalid origin',cprefix); %#ok<AGROW>
    end
    if isfield(metadata,'interpolation') && ~any(strcmp(char(string(metadata.interpolation)), ...
            {'linear','previous','none'}))
        errors{end+1} = sprintf('%s has invalid interpolation',cprefix); %#ok<AGROW>
    end
    if ~isfield(ch,'values') || ~isnumeric(ch.values) || ~isvector(ch.values)
        errors{end+1} = [cprefix '.values must be a numeric vector']; %#ok<AGROW>
        continue
    end
    if ~isfield(ch,'valid') || ~(islogical(ch.valid) || isnumeric(ch.valid)) || ...
            ~isvector(ch.valid) || numel(ch.valid) ~= numel(ch.values)
        errors{end+1} = [cprefix '.valid must be a vector matching values']; %#ok<AGROW>
        continue
    end
    valid = logical(ch.valid(:));
    values = double(ch.values(:));
    if any(~valid & ~isnan(values))
        errors{end+1} = sprintf('%s invalid samples must be NaN',cprefix); %#ok<AGROW>
    end
    if any(valid & ~isfinite(values))
        errors{end+1} = sprintf('%s valid samples must be finite',cprefix); %#ok<AGROW>
    end
    if isfield(metadata,'axis_id') && ~any(axisIds == string(metadata.axis_id))
        errors{end+1} = sprintf('%s refers to a missing axis',cprefix);
    elseif isfield(metadata,'axis_id')
        a = lap.axes(axisIds == string(metadata.axis_id));
        if numel(a.time_s) ~= numel(values)
            errors{end+1} = sprintf('%s length does not match axis',cprefix); %#ok<AGROW>
        end
    end
end
end

function [errors,n,ok] = checkFiniteMonotonicVector(s,prefix,name,errors,monotonic)
n = [];
ok = false;
if ~isfield(s,name) || ~isnumeric(s.(name)) || ~isvector(s.(name)) || isempty(s.(name))
    errors{end+1} = sprintf('%s.%s must be a non-empty numeric vector',prefix,name); %#ok<AGROW>
    return
end
v = double(s.(name)(:));
n = numel(v);
if any(~isfinite(v))
    errors{end+1} = sprintf('%s.%s must be finite',prefix,name); %#ok<AGROW>
elseif monotonic && any(diff(v) < 0)
    errors{end+1} = sprintf('%s.%s must be monotonic non-decreasing',prefix,name); %#ok<AGROW>
else
    ok = true;
end
end

function [ok,value] = textField(s,name)
ok = false;
value = "";
if ~isstruct(s) || ~isfield(s,name) || isempty(s.(name))
    return
end
candidate = string(s.(name));
if isscalar(candidate)
    value = candidate;
    ok = true;
end
end
