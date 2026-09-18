function result = buildTelemetryResult(source,opts)
%BUILDTELEMETRYRESULT Capture solved QSS outputs into the v1 MATLAB model.
%   SOURCE is read only; no event solver is called. It may be an existing
%   Events2 object, a canonical telemetry result, or the plain struct described
%   in QSS_TELEMETRY_USAGE.md. Native vectors are copied in full. Vectors with
%   different lengths receive separate axes instead of being trimmed.

if nargin < 2 || isempty(opts), opts = struct(); end
opts = telemetryOptions(opts);
source = unwrapSource(source);
if isCanonicalResult(source)
    result = normalizeCanonicalResult(source,opts);
elseif isEventsResult(source)
    result = resultFromEvents2(source,opts);
else
    result = resultFromPlainSource(source,opts);
end

if opts.reconstruct_detail
    runtimeCar = opts.car;
    if isempty(runtimeCar) && isstruct(source) && isfield(source,'car'), runtimeCar = source.car; end
    if isempty(runtimeCar) && isEventsResult(source), runtimeCar = source.car; end
    for ci = 1:numel(result.cases)
        caseCar = runtimeCar;
        if isempty(caseCar) && isfield(result.cases(ci),'runtime_car'), caseCar = result.cases(ci).runtime_car; end
        for li = 1:numel(result.cases(ci).laps)
            lap = result.cases(ci).laps(li);
            reconLap = lap;
            if ~isfield(reconLap,'curvature_per_m')
                ti = find(strcmp({result.tracks.id},char(string(lap.track_id))),1);
                if ~isempty(ti), reconLap.curvature_per_m = result.tracks(ti).curvature_per_m; end
            end
            detail = reconstructLapTelemetry(caseCar,reconLap,opts);
            result.cases(ci).laps(li) = appendDetail(lap,detail);
        end
    end
end
result.manifest = completeManifest(result.manifest,opts,result);
[ok,report] = validateQSSResults(result);
if ~ok
    error('buildTelemetryResult:invalidResult','Built result failed validation: %s',strjoin(report.errors,'; '));
end
end

function opts = telemetryOptions(opts)
if ~isstruct(opts) || numel(opts) ~= 1, error('buildTelemetryResult:badOptions','opts must be a scalar struct.'); end
if isfield(opts,'reconstruct') && ~isfield(opts,'reconstruct_detail'), opts.reconstruct_detail = opts.reconstruct; end
if isfield(opts,'detail_points') && ~isfield(opts,'reconstruct_detail'), opts.reconstruct_detail = opts.detail_points; end
defaults = struct('reconstruct_detail',false,'spacing_m',2,'min_speed_mps',1,'max_evaluations',1500, ...
    'case_id','baseline','case_label','Baseline','source_type','matlab_qss','run_uuid','', ...
    'producer','QSS MATLAB telemetry exporter','overwrite',false,'unique_name',true,'car',[]);
names = fieldnames(defaults);
for i = 1:numel(names)
    if ~isfield(opts,names{i}) || isempty(opts.(names{i})), opts.(names{i}) = defaults.(names{i}); end
end
opts.reconstruct_detail = logical(opts.reconstruct_detail);
validateattributes(opts.spacing_m,{'numeric'},{'scalar','real','finite','positive'},mfilename,'opts.spacing_m');
validateattributes(opts.min_speed_mps,{'numeric'},{'scalar','real','finite','nonnegative'},mfilename,'opts.min_speed_mps');
validateattributes(opts.max_evaluations,{'numeric'},{'scalar','real','finite','positive'},mfilename,'opts.max_evaluations');
opts.max_evaluations = min(floor(opts.max_evaluations),1500);
end

function source = unwrapSource(source)
for i = 1:3
    if isstruct(source) && isfield(source,'result') && isstruct(source.result) && isscalar(source.result)
        source = source.result;
    elseif isstruct(source) && isfield(source,'results') && isstruct(source.results) && isscalar(source.results) && ~isfield(source,'cases')
        source = source.results;
    else
        break
    end
end
end

function tf = isCanonicalResult(source)
tf = isstruct(source) && isscalar(source) && isfield(source,'tracks') && isfield(source,'cases') && ...
    (isfield(source,'manifest') || isfield(source,'format'));
end

function tf = isEventsResult(source)
tf = isobject(source) && isprop(source,'autocross') && isprop(source,'endurance');
end

function result = normalizeCanonicalResult(source,opts)
result = struct('manifest',struct(),'tracks',repmat(emptyTrack(),1,0), ...
    'cases',repmat(emptyCase(),1,0),'runtime_car',[]);
if isfield(source,'manifest') && isstruct(source.manifest), result.manifest = safeValue(source.manifest); end
if ~isfield(source,'tracks') || isempty(source.tracks), error('buildTelemetryResult:noTracks','A result needs a track.'); end
items = structItems(source.tracks,'track');
result.tracks = repmat(emptyTrack(),1,numel(items));
for i = 1:numel(items), result.tracks(i) = normalizeTrack(items{i},sprintf('track_%d',i)); end
items = structItems(source.cases,'case');
result.cases = repmat(emptyCase(),1,numel(items));
for i = 1:numel(items), result.cases(i) = normalizeCase(items{i},opts,sprintf('case_%d',i),result.tracks); end
if isempty(result.cases), error('buildTelemetryResult:noCases','A result needs a case.'); end
end

function result = resultFromEvents2(ev,opts)
result = struct('manifest',struct(),'tracks',repmat(emptyTrack(),1,0), ...
    'cases',repmat(emptyCase(),1,0),'runtime_car',[]);
trackItems = {};
if isprop(ev,'autocross_track') && ~isempty(ev.autocross_track)
    trackItems{end+1} = makeTrackSource('autocross',ev.autocross_track, ...
        struct('name','Autocross','geometry_source','curvature')); %#ok<AGROW>
end
if isprop(ev,'endurance_track') && ~isempty(ev.endurance_track)
    trackItems{end+1} = makeTrackSource('endurance',ev.endurance_track, ...
        struct('name','Endurance','geometry_source','curvature')); %#ok<AGROW>
end
if isprop(ev,'accel') && ~isempty(ev.accel)
    at = fieldVector(ev.accel,{'time_vec','time_s','time'});
    if ~isempty(at)
        lengthValue = optionOr(ev.eventParams,'accel_length',max(numel(at)-1,1));
        trackItems{end+1} = struct('id','accel','distance_m',linspace(0,lengthValue,numel(at)).', ...
            'curvature_per_m',zeros(numel(at),1),'metadata',struct('name','Acceleration','geometry_source','none')); %#ok<AGROW>
    end
end
result.tracks = repmat(emptyTrack(),1,numel(trackItems));
for i = 1:numel(trackItems), result.tracks(i) = normalizeTrack(trackItems{i},sprintf('track_%d',i)); end
caseSource = struct('id',char(string(optionOr(opts,'case_id','baseline'))), ...
    'metadata',struct('label',char(string(optionOr(opts,'case_label','Baseline')))), ...
    'setup',struct(),'laps',repmat(emptyPlainLap(),1,0),'runtime_car',[]);
if isprop(ev,'car') && ~isempty(ev.car), caseSource.setup = carSetup(ev.car); end
lapItems = {};
if isprop(ev,'autocross') && isstruct(ev.autocross) && ~isempty(ev.autocross)
    lapItems{end+1} = eventLapSource(ev.autocross,'autocross','autocross','flying',ev); %#ok<AGROW>
end
if isprop(ev,'endurance') && isstruct(ev.endurance) && ~isempty(ev.endurance)
    lapItems{end+1} = eventLapSource(ev.endurance,'endurance','endurance','representative',ev); %#ok<AGROW>
end
if isprop(ev,'accel') && isstruct(ev.accel) && ~isempty(ev.accel)
    lapItems{end+1} = eventLapSource(ev.accel,'accel','accel','straight',ev); %#ok<AGROW>
end
if isempty(lapItems), error('buildTelemetryResult:noSolvedEvents','Events2 contains no solved event outputs.'); end
caseSource.laps = [lapItems{:}];
result.cases(1) = normalizeCase(caseSource,opts,'baseline',result.tracks);
result.manifest = struct('source_type','Events2');
result.runtime_car = ev.car;
end

function p = eventLapSource(eventOutput,eventName,trackId,role,ev)
p = struct('id',eventName,'event',eventName,'role',role,'track_id',trackId, ...
    'metadata',struct('event',eventName,'role',role),'time_s',[], ...
    'distance_m',[],'lap_time_s',[],'channels',repmat(emptyPlainChannel(),1,0));
p.time_s = fieldVector(eventOutput,{'time_vec','time_s','time'});
if isempty(p.time_s) && isprop(ev,'times') && isstruct(ev.times) && isfield(ev.times,eventName)
    p.lap_time_s = ev.times.(eventName);
end

if isfield(eventOutput,'long_vel') && ~isempty(eventOutput.long_vel)
    p.channels(end+1) = plainChannel('speed_mps','Vehicle speed','m/s',eventOutput.long_vel,'simulation_output','linear');
elseif isfield(eventOutput,'long_vel_vector') && ~isempty(eventOutput.long_vel_vector)
    p.channels(end+1) = plainChannel('speed_mps','Vehicle speed','m/s',eventOutput.long_vel_vector,'simulation_output','linear');
end
if isfield(eventOutput,'long_accel') && ~isempty(eventOutput.long_accel)
    p.channels(end+1) = plainChannel('long_accel_mps2','Longitudinal acceleration','m/s^2',eventOutput.long_accel,'simulation_output','linear');
elseif isfield(eventOutput,'long_accel_vector') && ~isempty(eventOutput.long_accel_vector)
    p.channels(end+1) = plainChannel('long_accel_mps2','Longitudinal acceleration','m/s^2',eventOutput.long_accel_vector,'simulation_output','linear');
end
if isfield(eventOutput,'lat_accel') && ~isempty(eventOutput.lat_accel)
    p.channels(end+1) = plainChannel('lat_accel_mps2','Lateral acceleration','m/s^2',eventOutput.lat_accel,'simulation_output','linear');
elseif isfield(eventOutput,'lat_accel_vector') && ~isempty(eventOutput.lat_accel_vector)
    p.channels(end+1) = plainChannel('lat_accel_mps2','Lateral acceleration','m/s^2',eventOutput.lat_accel_vector,'simulation_output','linear');
end
if isfield(eventOutput,'num_upshifts'), p.metadata.num_upshifts = eventOutput.num_upshifts; end
if strcmp(eventName,'endurance')
    p.metadata.representative_lap = true;
    if isprop(ev,'eventParams') && isstruct(ev.eventParams) && isfield(ev.eventParams,'endurance_laps'), p.metadata.lap_count = ev.eventParams.endurance_laps; end
    if isfield(eventOutput,'lap_time'), p.metadata.raw_lap_time_s = eventOutput.lap_time; end
end
if isprop(ev,'times') && isstruct(ev.times) && isfield(ev.times,eventName), p.metadata.adjusted_event_time_s = ev.times.(eventName); end
end

function result = resultFromPlainSource(source,opts)
if ~isstruct(source) || ~isscalar(source), error('buildTelemetryResult:badSource','SOURCE must be a scalar struct or Events2 object.'); end
trackItems = {};
if isfield(source,'track') && ~isempty(source.track), trackItems = structItems(source.track,'track');
elseif isfield(source,'tracks') && ~isempty(source.tracks), trackItems = structItems(source.tracks,'track'); end
if isfield(source,'laps') && ~isempty(source.laps), lapItems = structItems(source.laps,'lap');
elseif isfield(source,'case') && isstruct(source.case) && isfield(source.case,'laps'), lapItems = structItems(source.case.laps,'lap');
else, lapItems = {}; end
if isempty(lapItems) && isfield(source,'autocross'), lapItems{1} = eventLapSource(source.autocross,'autocross','autocross','flying',source); end
if isempty(lapItems), error('buildTelemetryResult:noLaps','SOURCE contains no laps.'); end
if isempty(trackItems)
    for i = 1:numel(lapItems)
        d = fieldVector(lapItems{i},{'distance_m','distance','s'});
        if isempty(d), d = (0:max(numel(fieldVector(lapItems{i},{'time_s','time_vec'}))-1,1)).'; end
        k = fieldVector(lapItems{i},{'curvature_per_m','curvature','kappa_track'});
        if isempty(k), k = zeros(size(d)); end
        trackItems{end+1} = struct('id',fieldText(lapItems{i},'track_id',sprintf('track_%d',i)), ...
            'distance_m',d,'curvature_per_m',coordinateToLength(k,numel(d)), ...
            'metadata',struct('geometry_source','curvature')); %#ok<AGROW>
    end
end
result = struct('manifest',struct(),'tracks',repmat(emptyTrack(),1,numel(trackItems)), ...
    'cases',repmat(emptyCase(),1,0),'runtime_car',[]);
for i = 1:numel(trackItems), result.tracks(i) = normalizeTrack(trackItems{i},sprintf('track_%d',i)); end
caseSource = source;
caseSource.id = fieldText(source,'case_id',fieldText(source,'id',char(string(opts.case_id))));
if ~isfield(caseSource,'metadata'), caseSource.metadata = fieldOr(source,'case_metadata',struct('label',opts.case_label)); end
if ~isfield(caseSource,'setup'), caseSource.setup = fieldOr(source,'setup',struct()); end
caseSource.laps = [lapItems{:}];
result.cases(1) = normalizeCase(caseSource,opts,caseSource.id,result.tracks);
if isfield(source,'source_type'), result.manifest.source_type = source.source_type; end
if isfield(source,'car'), result.runtime_car = source.car; end
end

function c = normalizeCase(source,opts,defaultId,tracks)
c = emptyCase();
c.id = fieldText(source,'id',defaultId);
c.metadata = safeValue(fieldOr(source,'metadata',fieldOr(source,'case_metadata',struct('label',opts.case_label))));
c.setup = safeValue(fieldOr(source,'setup',struct()));
lapItems = {};
if isfield(source,'laps') && ~isempty(source.laps), lapItems = structItems(source.laps,'lap'); end
if isempty(lapItems) && isfield(source,'lap') && ~isempty(source.lap), lapItems = structItems(source.lap,'lap'); end
if isempty(lapItems), error('buildTelemetryResult:noLaps','Case %s contains no laps.',c.id); end
c.laps = repmat(emptyLap(),1,numel(lapItems));
for i = 1:numel(lapItems)
    li = lapItems{i};
    trackId = fieldText(li,'track_id','');
    if isempty(trackId), trackId = tracks(min(i,numel(tracks))).id; end
    idx = find(strcmp({tracks.id},trackId),1);
    if isempty(idx), idx = 1; trackId = tracks(idx).id; end
    c.laps(i) = normalizeLap(li,tracks(idx),trackId,sprintf('lap_%d',i));
end
if isfield(source,'runtime_car'), c.runtime_car = source.runtime_car; end
end

function lap = normalizeLap(source,track,trackId,defaultId)
if isfield(source,'axes') && isfield(source,'channels') && ~isempty(source.axes)
    lap = normalizeCanonicalLap(source,trackId,defaultId);
else
    lap = makeLap(source,track,trackId,defaultId);
end
end

function lap = normalizeCanonicalLap(source,trackId,defaultId)
lap = emptyLap();
lap.id = fieldText(source,'id',defaultId);
lap.track_id = fieldText(source,'track_id',trackId);
lap.metadata = safeValue(fieldOr(source,'metadata',struct()));
% The on-disk schema carries the association in lap metadata.  Keep it when
% normalising an already canonical result too, so Python can load files with
% more than one track without relying on a single-track fallback.
lap.metadata.track_id = lap.track_id;
axisItems = structItems(source.axes,'axis');
lap.axes = repmat(emptyAxis(),1,numel(axisItems));
for i = 1:numel(axisItems)
    a = axisItems{i};
    lap.axes(i) = struct('id',fieldText(a,'id',sprintf('axis_%d',i)), ...
        'time_s',double(fieldVector(a,{'time_s'})),'distance_m',double(fieldVector(a,{'distance_m'})));
end

channelItems = structItems(source.channels,'channel');
lap.channels = repmat(emptyChannel(),1,numel(channelItems));
for i = 1:numel(channelItems)
    lap.channels(i) = normalizeChannel(channelItems{i},lap.axes,sprintf('channel_%d',i));
end
if isfield(source,'diagnostics'), lap.diagnostics = safeValue(source.diagnostics); end
end

function lap = makeLap(source,track,trackId,defaultId)
lap = emptyLap();
lap.id = fieldText(source,'id',defaultId);
lap.track_id = trackId;
lap.metadata = safeValue(fieldOr(source,'metadata',struct()));
lap.metadata.event = fieldText(source,'event',fieldText(lap.metadata,'event',''));
lap.metadata.role = fieldText(source,'role',fieldText(lap.metadata,'role',''));
lap.metadata.track_id = trackId;
if isfield(source,'lap_time_s'), lap.metadata.raw_lap_time_s = source.lap_time_s; end
rawTime = fieldVector(source,{'time_s','time_vec','time','t'});
rawDistance = fieldVector(source,{'distance_m','distance','s','arclength'});
if isempty(rawDistance), rawDistance = track.distance_m; end
if isempty(rawTime)
    duration = scalarField(source,{'lap_time_s','duration_s','time_total'},[]);
    if isempty(duration), duration = max(numel(rawDistance)-1,0); end
    rawTime = linspace(0,duration,max(numel(rawDistance),1)).';
end
[rawTime,lap.metadata] = normalizeTime(rawTime,lap.metadata);
lap.metadata.raw_time_source_field = firstField(source,{'time_s','time_vec','time','t'},'derived');
lap.metadata.time_semantics = 'raw_native_time';
if ~isfield(lap.metadata,'adjusted_time_semantics'), lap.metadata.adjusted_time_semantics = 'none'; end
channelItems = plainChannels(source);
if isempty(channelItems), error('buildTelemetryResult:noChannels','Lap %s contains no numeric channels.',lap.id); end
lap.axes = repmat(emptyAxis(),1,0);
lap.channels = repmat(emptyChannel(),1,0);
mapping = repmat(emptyMapping(),1,0);
for i = 1:numel(channelItems)
    item = channelItems(i); values = double(item.values(:));
    if isempty(values), continue; end
    axisId = item.axis_id;
    if isempty(axisId), axisId = chooseAxisId(lap.axes,rawTime,rawDistance,numel(values)); end
    axisIdx = find(strcmp({lap.axes.id},axisId),1);
    if isempty(axisIdx)
        lap.axes(end+1) = struct('id',axisId,'time_s',coordinateToLength(rawTime,numel(values)), ...
            'distance_m',coordinateToLength(rawDistance,numel(values))); %#ok<AGROW>
    end
    lap.channels(end+1) = makeChannel(item,axisId); %#ok<AGROW>
    mapItem = emptyMapping();
    mapItem.channel_id = item.id; mapItem.source_field = item.source_field;
    mapItem.source_length = numel(values); mapItem.axis_id = axisId;
    mapItem.node_index = (1:numel(values)).';
    mapItem.segment_index = max((1:numel(values)).'-1,1);
    mapping(end+1) = mapItem; %#ok<AGROW>
end
if isempty(lap.axes)
    n = max(numel(rawTime),numel(rawDistance));
    lap.axes = struct('id','native','time_s',coordinateToLength(rawTime,n), ...
        'distance_m',coordinateToLength(rawDistance,n));
end
lap.metadata.node_segment_mapping = mapping;
if isfield(source,'diagnostics'), lap.diagnostics = safeValue(source.diagnostics); end
end

function items = plainChannels(source)
items = repmat(emptyPlainChannel(),1,0);
if isfield(source,'channels') && ~isempty(source.channels)
    raw = structItems(source.channels,'channel');
    for i = 1:numel(raw)
        item = raw{i};
        if ~isfield(item,'values'), continue; end
        p = emptyPlainChannel();
        p.id = fieldText(item,'id',sprintf('channel_%d',i)); p.label = fieldText(item,'label',p.id);
        p.unit = fieldText(item,'unit',''); p.values = item.values; p.valid = fieldOr(item,'valid',[]);
        p.axis_id = fieldText(item,'axis_id',''); p.origin = fieldText(item,'origin','simulation_output');
        p.interpolation = fieldText(item,'interpolation','linear'); p.description = fieldText(item,'description','');
        p.coordinate_frame = fieldText(item,'coordinate_frame','vehicle'); p.sign_convention = fieldText(item,'sign_convention','');
        p.source_field = fieldText(item,'source_field',p.id); items(end+1) = p; %#ok<AGROW>
    end
    return
end
catalog = { ...
    'speed_mps',{'speed_mps','long_vel','longitudinal_velocity','vCar'},'Vehicle speed','m/s','linear'; ...
    'long_accel_mps2',{'long_accel_mps2','long_accel','longitudinal_acceleration'},'Longitudinal acceleration','m/s^2','linear'; ...
    'lat_accel_mps2',{'lat_accel_mps2','lat_accel','lateral_acceleration'},'Lateral acceleration','m/s^2','linear'; ...
    'steering_angle_rad',{'steering_angle_rad','steer_angle','steering'},'Steering angle','rad','linear'; ...
    'control_demand',{'control_demand','throttle','signed_control_demand'},'Signed control demand','1','linear'; ...
    'gear',{'gear','current_gear'},'Gear','1','previous'; ...
    'curvature_per_m',{'curvature_per_m','curvature'},'Track curvature','1/m','linear'};
for i = 1:size(catalog,1)
    id = catalog{i,1}; names = catalog{i,2}; name = firstField(source,names,'');
    if isempty(name) || ~isfield(source,name) || ~isnumeric(source.(name)), continue; end
    p = emptyPlainChannel(); p.id = id; p.label = catalog{i,3}; p.unit = catalog{i,4};
    p.values = source.(name); p.origin = 'simulation_output'; p.interpolation = catalog{i,5}; p.source_field = name;
    if strcmp(id,'steering_angle_rad') && contains(name,'steer'), p.values = deg2rad(double(p.values)); end
    items(end+1) = p; %#ok<AGROW>
end
fields = fieldnames(source);
for i = 1:numel(fields)
    name = fields{i}; isKnown = any(cellfun(@(x) any(strcmp(name,x{1})),catalog(:,2)));
    excluded = any(strcmp(name,{'id','event','role','track_id','metadata','time_s','time_vec','distance_m', ...
        'distance','s','arclength','channels','lap_time_s','duration_s','time_total'}));
    if isKnown || excluded, continue; end
    if isnumeric(source.(name)) && isvector(source.(name)) && ~isempty(source.(name))
        p = emptyPlainChannel(); p.id = matlab.lang.makeValidName(name); p.label = name;
        p.values = source.(name); p.source_field = name; items(end+1) = p; %#ok<AGROW>
    end
end
end
function ch = normalizeChannel(source,axes,defaultId)
if isfield(source,'metadata') && isstruct(source.metadata), metadata = source.metadata; else, metadata = source; end
ch = emptyChannel();
ch.id = fieldText(metadata,'id',fieldText(source,'id',defaultId)); ch.label = fieldText(metadata,'label',ch.id);
ch.unit = fieldText(metadata,'unit',''); ch.axis_id = fieldText(metadata,'axis_id','');
ch.origin = fieldText(metadata,'origin','simulation_output'); ch.interpolation = fieldText(metadata,'interpolation','linear');
ch.description = fieldText(metadata,'description',''); ch.coordinate_frame = fieldText(metadata,'coordinate_frame','vehicle');
ch.sign_convention = fieldText(metadata,'sign_convention',''); ch.values = double(fieldVector(source,{'values'}));
ch.valid = logical(fieldOr(source,'valid',isfinite(ch.values)));
if isempty(ch.axis_id)
    axisIdx = find(arrayfun(@(a) numel(a.time_s)==numel(ch.values),axes),1);
    if isempty(axisIdx), axisIdx = 1; end
    ch.axis_id = axes(axisIdx).id;
end
ch.values(~ch.valid(:)) = NaN; ch.valid = ch.valid(:);
end

function ch = makeChannel(item,axisId)
ch = emptyChannel(); ch.id = item.id; ch.label = item.label; ch.unit = item.unit; ch.axis_id = axisId;
ch.origin = item.origin; ch.interpolation = item.interpolation; ch.description = item.description;
ch.coordinate_frame = item.coordinate_frame; ch.sign_convention = item.sign_convention; ch.values = double(item.values(:));
if isempty(item.valid), ch.valid = isfinite(ch.values); else, ch.valid = logical(item.valid(:)); end
ch.valid = ch.valid(:) & isfinite(ch.values); ch.values(~ch.valid) = NaN;
end

function lap = appendDetail(lap,detail)
if isempty(lap.axes), lap.axes = detail.axes; else, lap.axes(end+1) = detail.axes(1); end
if isempty(lap.channels), lap.channels = detail.channels;
else
    for i = 1:numel(detail.channels)
        idx = find(strcmp({lap.channels.id},detail.channels(i).id),1);
        if isempty(idx), lap.channels(end+1) = detail.channels(i); else, lap.channels(idx) = detail.channels(i); end
    end
end
lap.reconstruction = detail.diagnostics; lap.diagnostics = detail.diagnostics;
if ~isstruct(lap.metadata), lap.metadata = struct(); end
lap.metadata.reconstruction = detail.metadata;
end

function track = normalizeTrack(source,defaultId)
track = emptyTrack(); track.id = fieldText(source,'id',defaultId);
track.metadata = safeValue(fieldOr(source,'metadata',struct()));
track.distance_m = double(fieldVector(source,{'distance_m','distance','s','arclength'}));
track.curvature_per_m = double(fieldVector(source,{'curvature_per_m','curvature','kappa_track'}));
if isempty(track.distance_m), error('buildTelemetryResult:badTrack','Track %s has no distance.',track.id); end
if isempty(track.curvature_per_m), track.curvature_per_m = zeros(size(track.distance_m)); end
track.distance_m = track.distance_m(:); track.curvature_per_m = coordinateToLength(track.curvature_per_m,numel(track.distance_m));
if isfield(source,'x_m') && ~isempty(source.x_m), track.x_m = double(source.x_m(:)); end
if isfield(source,'y_m') && ~isempty(source.y_m), track.y_m = double(source.y_m(:)); end
end

function t = makeTrackSource(id,trackArray,metadata)
t = struct('id',id,'metadata',metadata,'distance_m',trackArray(1,:).','curvature_per_m',trackArray(2,:).');
end

function p = plainChannel(id,label,unit,values,origin,interpolation)
p = emptyPlainChannel(); p.id = id; p.label = label; p.unit = unit; p.values = values;
p.origin = origin; p.interpolation = interpolation; p.source_field = id;
end

function m = completeManifest(m,opts,result)
if ~isstruct(m), m = struct(); end
m = safeValue(m); m.format = 'qss-telemetry'; m.schema_version = '1.0';
if ~isfield(m,'run_uuid') || isempty(m.run_uuid), m.run_uuid = makeUuid(opts.run_uuid); end
if ~isfield(m,'creation_utc') || isempty(m.creation_utc), m.creation_utc = utcNow(); end
if ~isfield(m,'source_type') || isempty(m.source_type), m.source_type = opts.source_type; end
if ~isfield(m,'producer') || isempty(m.producer), m.producer = opts.producer; end
if ~isfield(m,'matlab_version') || isempty(m.matlab_version), m.matlab_version = version; end
if ~isfield(m,'git_revision') || isempty(m.git_revision), [m.git_revision,m.git_dirty] = gitInfo(); end
if ~isfield(m,'git_dirty') || isempty(m.git_dirty), [~,m.git_dirty] = gitInfo(); end
m.tracks = {result.tracks.id}; m.cases = {result.cases.id};
m.track_inventory = arrayfun(@(t) struct('id',t.id,'metadata',t.metadata),result.tracks,'UniformOutput',false);
m.case_inventory = arrayfun(@(c) struct('id',c.id,'metadata',c.metadata,'lap_ids',{{c.laps.id}}),result.cases,'UniformOutput',false);
end

function [revision,dirty] = gitInfo()
revision = 'unknown'; dirty = false;
try
    [status,out] = system('git rev-parse --short HEAD');
    if status == 0 && ~isempty(strtrim(out)), revision = strtrim(out); end
    [status,out] = system('git status --porcelain --untracked-files=no');
    dirty = status == 0 && ~isempty(strtrim(out));
catch
end
end

function value = makeUuid(candidate)
if ~isempty(candidate), value = char(string(candidate)); return; end
try
    value = char(java.util.UUID.randomUUID());
catch
    value = sprintf('run-%s-%06d',datestr(now,'yyyymmddTHHMMSS'),randi(999999));
end
end

function value = utcNow()
try
    value = char(datetime('now','TimeZone','UTC','Format','yyyy-MM-dd''T''HH:mm:ss.SSS''Z'''));
catch
    value = [datestr(now,'yyyy-mm-ddTHH:MM:SS') 'Z'];
end
end

function output = safeValue(input)
if isstruct(input)
    output = repmat(struct(),size(input)); fields = fieldnames(input);
    for n = 1:numel(input)
        for i = 1:numel(fields)
            try, output(n).(fields{i}) = safeValue(input(n).(fields{i}));
            catch, output(n).(fields{i}) = class(input(n).(fields{i}));
            end
        end
    end
elseif iscell(input)
    output = cell(size(input)); for i = 1:numel(input), output{i} = safeValue(input{i}); end
elseif isobject(input)
    output = class(input);
elseif isa(input,'function_handle')
    output = func2str(input);
elseif isstring(input)
    output = cellstr(input);
elseif ischar(input) || islogical(input)
    output = input;
elseif isnumeric(input)
    output = double(input); output(~isfinite(output)) = NaN;
else
    output = char(string(input));
end
end

function items = structItems(value,kind)
if isempty(value), items = {}; return; end
if iscell(value), items = value; return; end
if ~isstruct(value), error('buildTelemetryResult:bad%s',kind,'Expected struct data.'); end
if isscalar(value)
    names = fieldnames(value);
    if ~ismember('id',names) && ~ismember('values',names) && all(cellfun(@(n) isstruct(value.(n)),names))
        items = cell(1,numel(names));
        for i = 1:numel(names)
            item = value.(names{i}); if ~isfield(item,'id'), item.id = names{i}; end
            items{i} = item;
        end
    else
        items = {value};
    end
else
    items = num2cell(value);
end
end

function value = fieldOr(s,name,default)
if isstruct(s) && isfield(s,name), value = s.(name); else, value = default; end
end

function value = optionOr(s,name,default)
if isstruct(s) && isfield(s,name), value = s.(name); else, value = default; end
end

function value = fieldText(s,name,default)
value = default;
if isstruct(s) && isfield(s,name) && ~isempty(s.(name)), value = char(string(s.(name))); end
end

function value = firstField(s,names,default)
value = default;
for i = 1:numel(names)
    if isstruct(s) && isfield(s,names{i}) && ~isempty(s.(names{i})), value = names{i}; return; end
end
end

function value = scalarField(s,names,default)
value = default;
for i = 1:numel(names)
    if isstruct(s) && isfield(s,names{i}) && isnumeric(s.(names{i})) && isscalar(s.(names{i})), value = double(s.(names{i})); return; end
end
end

function value = fieldVector(s,names)
value = [];
for i = 1:numel(names)
    if isstruct(s) && isfield(s,names{i}) && isnumeric(s.(names{i})) && isvector(s.(names{i}))
        value = double(s.(names{i})(:)); return
    end
end
end

function [value,metadata] = normalizeTime(value,metadata)
value = double(value(:)); if isempty(value), return; end
if any(~isfinite(value)), value = fillmissing(value,'linear','EndValues','nearest'); end
if any(diff(value)<0)
    if all(value>=0), value = cumsum(value); metadata.adjusted_time_semantics = 'cumulative_from_native_segments';
    else, value = cummax(value); metadata.adjusted_time_semantics = 'monotonic_clamped'; end
end
if value(1)<0, value = value-value(1); end
end

function value = coordinateToLength(coord,n)
coord = double(coord(:));
if n<=0, value=zeros(0,1);
elseif isempty(coord), value=(0:n-1).';
elseif numel(coord)==n, value=coord;
elseif n==1, value=coord(1);
else, value=linspace(coord(1),coord(end),n).'; end
if any(~isfinite(value)), value=fillmissing(value,'linear','EndValues','nearest'); end
if any(diff(value)<0), value=cummax(value); end
end

function id = chooseAxisId(axes,time,distance,n)
if isempty(axes), id='native'; return; end
mask = arrayfun(@(a) numel(a.time_s)==n && numel(a.distance_m)==n,axes);
if any(mask), id=axes(find(mask,1)).id; else, id=sprintf('native_%d',n); end
end

function setup = carSetup(car)
setup = struct('class',class(car)); fields={'M','W_b','l_f','l_r','t_f','t_r','h_g','R','I_zz','g','Crr'};
for i=1:numel(fields)
    try, if isprop(car,fields{i}), setup.(fields{i})=double(car.(fields{i})); end, catch, end
end
end

function c = emptyCase()
c = struct('id','','metadata',struct(),'setup',struct(),'laps',repmat(emptyLap(),1,0),'runtime_car',[]);
end

function t = emptyTrack()
t = struct('id','','metadata',struct(),'distance_m',zeros(0,1),'curvature_per_m',zeros(0,1),'x_m',[],'y_m',[]);
end

function l = emptyLap()
l = struct('id','','metadata',struct(),'track_id','','axes',repmat(emptyAxis(),1,0), ...
    'channels',repmat(emptyChannel(),1,0),'diagnostics',struct(),'reconstruction',struct());
end

function a = emptyAxis()
a = struct('id','','time_s',zeros(0,1),'distance_m',zeros(0,1));
end

function c = emptyChannel()
c = struct('id','','label','','unit','','axis_id','','origin','simulation_output','interpolation','linear', ...
    'description','','coordinate_frame','vehicle','sign_convention','','values',zeros(0,1),'valid',false(0,1));
end

function p = emptyPlainChannel()
p = struct('id','','label','','unit','','axis_id','','origin','simulation_output','interpolation','linear', ...
    'description','','coordinate_frame','vehicle','sign_convention','','values',zeros(0,1),'valid',[],'source_field','');
end

function p = emptyPlainLap()
p = struct('id','','event','','role','','track_id','','metadata',struct(),'time_s',[],'distance_m',[],'lap_time_s',[], ...
    'channels',repmat(emptyPlainChannel(),1,0));
end

function m = emptyMapping()
m = struct('channel_id','','source_field','','source_length',0,'axis_id','','node_index',zeros(0,1),'segment_index',zeros(0,1));
end
