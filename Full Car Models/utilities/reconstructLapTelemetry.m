function detail = reconstructLapTelemetry(car,lap,opts)
%RECONSTRUCTLAPTELEMETRY Reconstruct validated quasi-steady detail points.
%   The returned detail object has one qss_detail axis, qss state channels,
%   numeric diagnostic channels, and a diagnostics struct with status and
%   failure reasons. Failed points remain NaN/invalid; no interpolation is
%   performed. QSS state order is
%   [steering_rad, signed_control_demand, v_long, v_lat, yaw_rate,
%    kappa_FL, kappa_FR, kappa_RL, kappa_RR].
%
%   Existing lap arrays are only sampled at selected points. The solver is
%   never called by this function for native capture.

if nargin < 3 || isempty(opts), opts = struct(); end
opts = reconOptions(opts);
[distance,time,curvature,speed,longAccel,latAccel,gear,shiftMask] = lapInputs(lap);
indices = selectPoints(distance,curvature,gear,shiftMask,longAccel,opts.spacing_m);
n = numel(indices);
detail = struct();
detail.axes = struct('id','qss_detail','time_s',time(indices),'distance_m',distance(indices));
detail.channels = repmat(emptyChannel(),1,0);
detail.diagnostics = emptyDiagnostics(n);
detail.metadata = struct('state_order',{{'steering_rad','signed_control_demand', ...
    'longitudinal_velocity_mps','lateral_velocity_mps','yaw_rate_radps', ...
    'kappa_fl','kappa_fr','kappa_rl','kappa_rr'}}, ...
    'spacing_m',opts.spacing_m,'max_evaluations',opts.max_evaluations, ...
    'residual_limits',residualLimits(),'selected_source_indices',indices(:), ...
    'unsupported_states',{{'low_speed_launch','shift'}});

stateIds = {'qss_steering_rad','qss_control_demand','qss_longitudinal_velocity_mps', ...
    'qss_lateral_velocity_mps','qss_yaw_rate_radps','qss_kappa_fl','qss_kappa_fr', ...
    'qss_kappa_rl','qss_kappa_rr'};
stateUnits = {'rad','1','m/s','m/s','rad/s','1','1','1','1'};
stateLabels = {'QSS steering','QSS signed control demand','QSS longitudinal velocity', ...
    'QSS lateral velocity','QSS yaw rate','QSS FL slip ratio','QSS FR slip ratio', ...
    'QSS RL slip ratio','QSS RR slip ratio'};
for i = 1:numel(stateIds)
    detail.channels(end+1) = makeDetailChannel(stateIds{i},stateLabels{i},stateUnits{i},detail.axes.id,n); %#ok<AGROW>
end

prior = [];
metricStore = struct();
for j = 1:n
    sourceIndex = indices(j);
    v = scalarAt(speed,sourceIndex,distance,detail.axes.distance_m(j));
    k = scalarAt(curvature,sourceIndex,distance,detail.axes.distance_m(j));
    ax = scalarAt(longAccel,sourceIndex,distance,detail.axes.distance_m(j));
    ay = scalarAt(latAccel,sourceIndex,distance,detail.axes.distance_m(j));
    if ~isfinite(ay) && isfinite(v) && isfinite(k), ay = v^2*k; end
    if ~isfinite(ax), ax = 0; end
    if ~isfinite(v), detail.diagnostics=setUnavailable(detail.diagnostics,j,'missing_velocity'); continue; end
    mode = 'propulsion';
    if ax < -1e-6, mode = 'braking'; end
    detail.diagnostics.mode{j} = mode;
    detail.diagnostics.distance_m(j) = detail.axes.distance_m(j);
    detail.diagnostics.time_s(j) = detail.axes.time_s(j);
    detail.diagnostics.target_long_accel_mps2(j) = ax;
    detail.diagnostics.target_lat_accel_mps2(j) = ay;
    if v <= opts.min_speed_mps && ax > 0
        detail.diagnostics=setUnavailable(detail.diagnostics,j,'low_speed_launch_unsupported'); continue;
    end
    if shiftMaskAt(shiftMask,sourceIndex,distance,detail.axes.distance_m(j))
        detail.diagnostics=setUnavailable(detail.diagnostics,j,'shift_unsupported'); continue;
    end
    if ~carAvailable(car)
        detail.diagnostics=setUnavailable(detail.diagnostics,j,'car_model_unavailable'); continue;
    end
    targetYaw = 0;
    if isfinite(v) && isfinite(k), targetYaw = v*k; end
    [P,solveInfo] = solvePoint(car,v,ay,ax,targetYaw,mode,prior,opts);
    detail.diagnostics.eval_count(j) = solveInfo.eval_count;
    detail.diagnostics.acceleration_residual_mps2(j) = solveInfo.accel_residual;
    detail.diagnostics.yaw_acceleration_residual_radps2(j) = solveInfo.yaw_residual;
    detail.diagnostics.torque_residual_Nm(j) = solveInfo.torque_residual;
    detail.diagnostics.min_virtual_load_N(j) = solveInfo.min_virtual_load;
    detail.diagnostics.aero_residual_in(j) = solveInfo.aero_residual;
    detail.diagnostics.failure_reason{j} = solveInfo.failure_reason;
    if ~solveInfo.valid
        detail.diagnostics.status{j} = 'failed';
        continue
    end
    detail.diagnostics.status{j} = 'solved';
    detail.diagnostics.valid(j) = true;
    prior = P;
    state = [deg2rad(P(1)),P(2),P(3),P(4),P(5),P(6:9)];
    for i = 1:numel(stateIds), detail.channels(i).values(j) = state(i); detail.channels(i).valid(j) = true; end
    metricStore = storeMetrics(metricStore,solveInfo.metrics,j,n);
end
detail.channels = addMetricChannels(detail.channels,metricStore,detail.axes.id,n);
detail.channels = addDiagnosticChannels(detail.channels,detail.diagnostics,detail.axes.id,n);
detail.metadata.status_codes = struct('unavailable',0,'solved',1,'failed',2);
detail.metadata.valid_points = nnz(detail.diagnostics.valid);
end

function opts = reconOptions(opts)
if ~isstruct(opts) || numel(opts)~=1, error('reconstructLapTelemetry:badOptions','opts must be a scalar struct.'); end
defaults = struct('spacing_m',2,'min_speed_mps',1,'max_evaluations',1500, ...
    'min_virtual_load_N',-0.1,'max_accel_residual',1e-3,'max_yaw_residual',1e-3, ...
    'max_torque_residual',1e-2,'max_aero_residual_in',1e-7,'solver_display','off');
names = fieldnames(defaults);
for i=1:numel(names), if ~isfield(opts,names{i}) || isempty(opts.(names{i})), opts.(names{i})=defaults.(names{i}); end, end
opts.max_evaluations = min(floor(opts.max_evaluations),1500);
validateattributes(opts.spacing_m,{'numeric'},{'scalar','real','finite','positive'},mfilename,'opts.spacing_m');
validateattributes(opts.min_speed_mps,{'numeric'},{'scalar','real','finite','nonnegative'},mfilename,'opts.min_speed_mps');
end

function limits = residualLimits()
limits = struct('acceleration_mps2',1e-3,'yaw_acceleration_radps2',1e-3, ...
    'torque_Nm',1e-2,'virtual_load_N',-0.1,'aero_residual_in',1e-7,'max_evaluations',1500);
end

function [distance,time,curvature,speed,longAccel,latAccel,gear,shiftMask] = lapInputs(lap)
[time,distance] = firstAxis(lap);
curvature = fieldVector(lap,{'curvature_per_m','curvature'});
speed = channelVector(lap,{'speed_mps','long_vel','vCar','qss_longitudinal_velocity_mps'});
longAccel = channelVector(lap,{'long_accel_mps2','long_accel','gLong'});
latAccel = channelVector(lap,{'lat_accel_mps2','lat_accel','gLat'});
gear = channelVector(lap,{'gear','current_gear'});
shiftMask = channelVector(lap,{'shift','shifting','shift_flag'});
if isempty(curvature), curvature = channelVector(lap,{'curvature_per_m'}); end
distance = distance(:); time = time(:);
if isempty(curvature), curvature = zeros(size(distance)); else, curvature = alignVector(curvature,numel(distance)); end
if isempty(speed), speed = nan(size(distance)); else, speed = alignVector(speed,numel(distance)); end
if isempty(longAccel), longAccel = nan(size(distance)); else, longAccel = alignVector(longAccel,numel(distance)); end
if isempty(latAccel), latAccel = nan(size(distance)); else, latAccel = alignVector(latAccel,numel(distance)); end
if isempty(gear), gear = nan(size(distance)); else, gear = alignVector(gear,numel(distance)); end
if isempty(shiftMask), shiftMask = false(size(distance)); else, shiftMask = alignVector(shiftMask,numel(distance))~=0; end
end

function [time,distance] = firstAxis(lap)
time = []; distance = [];
if isfield(lap,'axes') && ~isempty(lap.axes)
    time = double(lap.axes(1).time_s(:)); distance = double(lap.axes(1).distance_m(:));
end
if isempty(time), time = fieldVector(lap,{'time_s','time_vec','time','t'}); end
if isempty(distance), distance = fieldVector(lap,{'distance_m','distance','s','arclength'}); end
n = max(numel(time),numel(distance));
if n==0, n=1; end
if isempty(time), time=(0:n-1).'; else, time=alignVector(time,n); end
if isempty(distance), distance=(0:n-1).'; else, distance=alignVector(distance,n); end
time = monotonicCoordinate(time); distance = monotonicCoordinate(distance);
end

function values = channelVector(lap,ids)
values = [];
if isfield(lap,'channels') && isstruct(lap.channels)
    for i=1:numel(ids)
        idx=find(strcmp({lap.channels.id},ids{i}),1);
        if ~isempty(idx), values=double(lap.channels(idx).values(:)); return; end
    end
end
values = fieldVector(lap,ids);
if isempty(values) && isfield(lap,'metrics') && isstruct(lap.metrics)
    for i=1:numel(ids)
        if isfield(lap.metrics,ids{i}) && isnumeric(lap.metrics.(ids{i})), values=double(lap.metrics.(ids{i})(:)); return; end
    end
end
end

function values = fieldVector(s,names)
values = [];
for i=1:numel(names)
    if isstruct(s) && isfield(s,names{i}) && isnumeric(s.(names{i})) && isvector(s.(names{i}))
        values=double(s.(names{i})(:)); return
    end
end
end

function values = alignVector(values,n)
values = double(values(:));
if numel(values)==n, return; end
if isempty(values), values=nan(n,1); return; end
if n==1, values=values(1); return; end
if numel(values)==1, values=repmat(values,n,1); return; end
values=interp1(linspace(0,1,numel(values)),values,linspace(0,1,n),'nearest','extrap').';
end

function values = monotonicCoordinate(values)
values=double(values(:));
if isempty(values), values=0; end
if any(~isfinite(values)), values=fillmissing(values,'linear','EndValues','nearest'); end
if any(diff(values)<0), values=cummax(values); end
if values(1)<0, values=values-values(1); end
end

function indices = selectPoints(distance,curvature,gear,shiftMask,longAccel,spacing)
n=numel(distance); if n==0, indices=1; return; end
D=distance(end)-distance(1); targets=distance(1):spacing:distance(end); %#ok<NASGU>
indices=1;
for i=1:numel(targets), [~,k]=min(abs(distance-targets(i))); indices(end+1)=k; end %#ok<AGROW>
indices(end+1)=n;
if numel(curvature)>2
    slope=diff(curvature);
    for i=2:numel(curvature)-1
        if (slope(i-1)>=0 && slope(i)<=0) || (slope(i-1)<=0 && slope(i)>=0)
            indices(end+1:end+2)=[i i+1]; %#ok<AGROW>
        end
    end
end
transition = [find(diff(gear)~=0);find(diff(shiftMask)~=0);find(diff(signTransition(longAccel))~=0)];
for i=1:numel(transition)
    k=transition(i)+1; indices(end+1:end+2)=[max(k-1,1) min(k+1,n)]; %#ok<AGROW>
end
indices=unique(indices,'stable'); indices=indices(:);
end

function values = signTransition(values)
values=double(values(:)); values(~isfinite(values))=0; values=sign(values); values(abs(values)<1e-6)=0;
end

function value = scalarAt(values,index,coordinate,target)
value=NaN;
if isempty(values), return; end
if index<=numel(values), value=double(values(index)); end
if ~isfinite(value) && numel(values)>1 && ~isempty(coordinate)
    q=linspace(coordinate(1),coordinate(end),numel(values));
    value=interp1(q,values,target,'nearest','extrap');
end
end

function tf = shiftMaskAt(mask,index,coordinate,target)
tf=false;
if isempty(mask), return; end
if index<=numel(mask), tf=logical(mask(index)); end
if ~tf && numel(mask)>1 && index>1 && mask(index)~=mask(index-1), tf=true; end
if ~tf && numel(mask)>1 && ~isempty(coordinate)
    q=linspace(coordinate(1),coordinate(end),numel(mask)); tf=any(abs(q-target)<=max(diff(q),[],'omitnan')/2 & mask);
end
end

function tf = carAvailable(car)
tf=false;
if isobject(car), tf=ismethod(car,'equations') && ismethod(car,'metrics');
elseif isstruct(car), tf=isfield(car,'equations') && isa(car.equations,'function_handle') && isfield(car,'metrics') && isa(car.metrics,'function_handle'); end
end

function [P,info] = solvePoint(car,v,targetAy,targetAx,targetYaw,mode,prior,opts)
info=emptySolveInfo(); P=[]; info.failure_reason='no_feasible_candidate';
if ~isfinite(v) || v<=0, info.failure_reason='invalid_velocity'; return; end
[bounds,sourceBounds]=stateBounds(v,targetYaw,mode);
seeds=candidateSeeds(car,v,targetAy,targetAx,targetYaw,mode,prior,bounds);
best=emptyCandidate(); budget=opts.max_evaluations;
for i=1:numel(seeds)
    if budget<=0, break; end
    seed=clipState(seeds{i},bounds);
    candidate=emptyCandidate(); candidate.P=seed; candidate.eval_count=0;
    if exist('fmincon','file')==2
        try
            maxEval=max(25,min(ceil(budget),opts.max_evaluations));
            fopts=optimoptions('fmincon','Display',opts.solver_display,'MaxFunctionEvaluations',maxEval, ...
                'MaxIterations',maxEval,'ConstraintTolerance',1e-8,'OptimalityTolerance',1e-8,'StepTolerance',1e-10);
            objective=@(x) stateObjective(x,seed,prior);
            constraint=@(x) qssConstraint(car,x,targetAy,targetAx,targetYaw);
            [x,~,flag,out]=fmincon(objective,seed,[],[],[],[],bounds.lb,bounds.ub,constraint,fopts);
            candidate.P=x; candidate.eval_count=fieldOr(out,'funcCount',maxEval); candidate.exitflag=flag;
        catch err
            candidate.failure_reason=['fmincon:' err.identifier];
        end
    end
    ev=evaluateState(car,candidate.P,targetAy,targetAx,targetYaw);
    candidate=scoreCandidate(candidate,ev,bounds,opts);
    budget=budget-candidate.eval_count;
    if betterCandidate(candidate,best), best=candidate; end
end
if ~isempty(best.P)
    P=best.P; info=best.info; info.eval_count=sum([best.eval_count]); info.metrics=best.metrics;
    info.source_bounds=sourceBounds;
end
end

function [bounds,sourceBounds] = stateBounds(v,targetYaw,mode)
if strcmp(mode,'braking')
    steerLimit=22; throttle=[-1,0]; kappa=[-0.2,0;-0.2,0;-0.2,0;-0.2,0];
else
    steerLimit=25; throttle=[0,1]; kappa=[0,0;0,0;0,0.2;0,0.2];
end
bounds.lb=[-steerLimit,throttle(1),v,-3,max(-2,targetYaw),kappa(:,1).'];
bounds.ub=[ steerLimit,throttle(2),v, 3,min(2,targetYaw),kappa(:,2).'];
sourceBounds=struct('steering_deg',[-steerLimit steerLimit],'control_demand',throttle, ...
    'longitudinal_velocity_mps',[v v],'lateral_velocity_mps',[-3 3], ...
    'yaw_rate_radps',[max(-2,targetYaw) min(2,targetYaw)],'kappa',kappa);
end

function seeds = candidateSeeds(car,v,targetAy,targetAx,targetYaw,mode,prior,bounds)
seeds={};
if ~isempty(prior), seeds{end+1}=prior; end
gg=nearbyGgSeed(car,v,targetAx,targetAy,mode);
if ~isempty(gg), seeds{end+1}=gg; end
seeds{end+1}=analyticalSeed(car,v,targetAx,targetYaw,mode);
seeds{end+1}=neutralSeed(v,targetYaw,mode);
for i=1:numel(seeds), seeds{i}=clipState(seeds{i},bounds); end
end

function seed = analyticalSeed(car,v,targetAx,targetYaw,mode)
steer=0; latVel=0; drag=0; mass=1; wheelRadius=0.3;
try, steer=atan2d(car.W_b*targetYaw/v,1); latVel=car.l_r*targetYaw; mass=car.M; wheelRadius=car.R; catch, end %#ok<NASGU>
try, drag=car.aero.drag(v); catch, end
demand=(targetAx+drag/max(mass,eps))/10;
if strcmp(mode,'braking'), demand=min(demand,0); else, demand=max(demand,0); end
if strcmp(mode,'braking'), k=[-0.002 -0.002 -0.002 -0.002]; else, k=[0 0 0.01 0.01]; end
seed=[steer,demand,v,latVel,targetYaw,k];
end

function seed = neutralSeed(v,targetYaw,mode)
if strcmp(mode,'braking'), demand=-0.01; k=[-0.001 -0.001 -0.001 -0.001]; else, demand=0.01; k=[0 0 0.005 0.005]; end
seed=[0,demand,v,0,targetYaw,k];
end

function seed = nearbyGgSeed(car,v,targetAx,targetAy,mode)
seed=[];
names={'ss_info','accel_info','decel_info'};
if strcmp(mode,'propulsion'), names={'accel_info','ss_info'}; elseif strcmp(mode,'braking'), names={'decel_info','ss_info'}; end
for n=1:numel(names)
    try, rows=car.(names{n}); catch, rows=[]; end
    if isempty(rows) || ~isnumeric(rows) || size(rows,2)<12, continue; end
    % x/P occupies columns 4:12 in the native g-g solution records.
    vcol=rows(:,6); [~,idx]=min(abs(vcol-v));
    seed=rows(idx,4:12);
    if numel(seed)==9, return; end
end
end

function state=clipState(state,bounds)
state=double(state(:).'); state=min(max(state,bounds.lb),bounds.ub);
end

function value=stateObjective(x,seed,prior)
weights=[1e-3 1 1e-2 1e-2 1e-2 ones(1,4)];
value=sum(weights.*(x(:).'-seed).^2);
if ~isempty(prior), value=value+1e-2*sum((x(:).'-prior(:).').^2); end
end

function [c,ceq]=qssConstraint(car,P,targetAy,targetAx,targetYaw)
ev=evaluateState(car,P,targetAy,targetAx,targetYaw);
c=-ev.Fzvirtual(:)-0.1;
ceq=[ev.lat_residual,ev.long_residual,ev.yaw_residual,ev.wheel_residual(:).'];
if any(~isfinite([c(:);ceq(:)])), c=1e12; ceq=1e12; end
end

function ev=evaluateState(car,P,targetAy,targetAx,targetYaw)
ev=struct('lat_residual',NaN,'long_residual',NaN,'yaw_residual',NaN,'wheel_residual',nan(4,1), ...
    'Fzvirtual',nan(4,1),'aero_residual',NaN,'metrics',struct());
try
    if isobject(car)
        [~,~,ev.lat_residual,actualAx,ev.yaw_residual,ev.wheel_residual,~,~,ev.Fzvirtual,~,~,~,~,~,~,ss]=car.equations(P);
        ev.long_residual=actualAx-targetAx;
        if isstruct(ss) && isfield(ss,'aero_residual_in'), ev.aero_residual=ss.aero_residual_in; end
        ev.metrics=car.metrics(P);
    else
        out=car.equations(P); if isstruct(out), ev=mergeEvaluation(ev,out,targetAx); end
        ev.metrics=car.metrics(P);
    end
catch
end
if isempty(ev.metrics), ev.metrics=struct(); end
end

function ev=mergeEvaluation(ev,out,targetAx)
names={'lat_residual','lat_accel','long_residual','long_accel','yaw_residual','yaw_accel','wheel_residual','wheel_accel','Fzvirtual','aero_residual'};
for i=1:numel(names)
    if isfield(out,names{i})
        switch names{i}
            case {'lat_residual','lat_accel'}, ev.lat_residual=out.(names{i});
            case {'long_residual','long_accel'}, ev.long_residual=out.(names{i})-targetAx;
            case {'yaw_residual','yaw_accel'}, ev.yaw_residual=out.(names{i});
            case {'wheel_residual','wheel_accel'}, ev.wheel_residual=out.(names{i});
            case 'Fzvirtual', ev.Fzvirtual=out.(names{i});
            case 'aero_residual', ev.aero_residual=out.(names{i});
        end
    end
end
end

function candidate=scoreCandidate(candidate,ev,bounds,opts)
candidate.metrics=ev.metrics; candidate.ev=ev;
candidate.info=emptySolveInfo();
candidate.info.metrics=ev.metrics;
candidate.info.accel_residual=max(abs([ev.lat_residual,ev.long_residual]));
candidate.info.yaw_residual=abs(ev.yaw_residual);
candidate.info.torque_residual=max(abs(ev.wheel_residual));
candidate.info.min_virtual_load=min(ev.Fzvirtual);
candidate.info.aero_residual=abs(ev.aero_residual);
candidate.info.failure_reason='residual_limits';
candidate.info.valid=all(isfinite([candidate.info.accel_residual,candidate.info.yaw_residual,candidate.info.torque_residual,candidate.info.min_virtual_load])) && ...
    candidate.info.accel_residual<=opts.max_accel_residual && candidate.info.yaw_residual<=opts.max_yaw_residual && ...
    candidate.info.torque_residual<=opts.max_torque_residual && candidate.info.min_virtual_load>=opts.min_virtual_load_N && ...
    (isnan(candidate.info.aero_residual) || candidate.info.aero_residual<=opts.max_aero_residual_in) && ...
    all(candidate.P>=bounds.lb-1e-10 & candidate.P<=bounds.ub+1e-10);
if ~candidate.info.valid
    if ~all(isfinite(candidate.P)), candidate.info.failure_reason='nonfinite_solution'; end
end
end

function tf=betterCandidate(a,b)
if isempty(b.P), tf=true; return; end
if a.info.valid~=b.info.valid, tf=a.info.valid; return; end
tf=maxFinite([a.info.accel_residual,a.info.yaw_residual,a.info.torque_residual]) < ...
    maxFinite([b.info.accel_residual,b.info.yaw_residual,b.info.torque_residual]);
end

function value=maxFinite(values)
values=values(isfinite(values)); if isempty(values), value=Inf; else, value=max(values); end
end

function info=emptySolveInfo()
info=struct('valid',false,'failure_reason','','eval_count',0,'accel_residual',NaN,'yaw_residual',NaN, ...
    'torque_residual',NaN,'min_virtual_load',NaN,'aero_residual',NaN,'metrics',struct(),'source_bounds',struct());
end

function candidate=emptyCandidate()
candidate=struct('P',[],'eval_count',0,'exitflag',0,'info',emptySolveInfo(),'metrics',struct(),'ev',struct());
end

function diagnostics=emptyDiagnostics(n)
% Assign fields after constructing the scalar struct.  Passing cell arrays
% directly to struct(...) expands them into a struct array, which prevents
% indexed updates such as diagnostics.mode{j} below.
diagnostics = struct();
diagnostics.status = repmat({'unavailable'},n,1);
diagnostics.valid = false(n,1);
diagnostics.mode = repmat({'unknown'},n,1);
diagnostics.failure_reason = repmat({'not_evaluated'},n,1);
diagnostics.eval_count = nan(n,1);
diagnostics.acceleration_residual_mps2 = nan(n,1);
diagnostics.yaw_acceleration_residual_radps2 = nan(n,1);
diagnostics.torque_residual_Nm = nan(n,1);
diagnostics.min_virtual_load_N = nan(n,1);
diagnostics.aero_residual_in = nan(n,1);
diagnostics.target_long_accel_mps2 = nan(n,1);
diagnostics.target_lat_accel_mps2 = nan(n,1);
diagnostics.distance_m = nan(n,1);
diagnostics.time_s = nan(n,1);
end

function diagnostics=setUnavailable(diagnostics,j,reason)
diagnostics.status{j}='unavailable'; diagnostics.valid(j)=false; diagnostics.failure_reason{j}=reason;
end

function ch=makeDetailChannel(id,label,unit,axisId,n)
ch=struct('id',id,'label',label,'unit',unit,'axis_id',axisId,'origin','qss_reconstructed', ...
    'interpolation','none','description','QSS detail channel','coordinate_frame','vehicle', ...
    'sign_convention','','values',nan(n,1),'valid',false(n,1));
end

function store=storeMetrics(store,metrics,index,n)
if ~isstruct(metrics) || isempty(metrics), return; end
names=fieldnames(metrics);
for i=1:numel(names)
    name=names{i}; value=metrics.(name);
    if isnumeric(value) && isscalar(value)
        if ~isfield(store,name), store.(name)=nan(n,1); end
        if isfinite(value), store.(name)(index)=double(value); end
    end
end
end

function channels=addMetricChannels(channels,store,axisId,n)
if isempty(fieldnames(store)), return; end
names=fieldnames(store);
for i=1:numel(names)
    name=names{i}; id=name;
    if any(strcmp({channels.id},id)), id=['metric_' name]; end
    [label,unit,interp,signConvention]=metricMetadata(name);
    ch=makeDetailChannel(id,label,unit,axisId,n); ch.values=store.(name); ch.valid=isfinite(ch.values); ch.description='Car.metrics field'; ch.interpolation=interp; ch.sign_convention=signConvention;
    channels(end+1)=ch; %#ok<AGROW>
end
end

function [label,unit,interp,signConvention]=metricMetadata(name)
label=name; unit=''; interp='linear'; signConvention='';
switch name
    case {'engine_rpm'}, label='Engine speed'; unit='rpm'; interp='linear';
    case {'current_gear'}, label='Gear'; unit='1'; interp='previous';
    case {'gLat','gLong','gMag'}, unit='g';
    case {'steer_angle'}, unit='deg'; signConvention='left positive';
    case {'throttle'}, unit='1'; signConvention='braking negative, propulsion positive';
    case {'vCar','lat_vel'}, unit='m/s';
    case {'yaw_rate'}, unit='rad/s'; signConvention='left positive';
    case {'downforce','drag','aero_downforce_front_N','aero_downforce_rear_N'}, unit='N';
    case {'aero_residual_in','front_ride_height_in','rear_ride_height_in','ride_height_FL_in','ride_height_FR_in','ride_height_RL_in','ride_height_RR_in'}, unit='in';
    case {'Fz_1','Fz_2','Fz_3','Fz_4','Fx_1','Fx_2','Fx_3','Fx_4','Fy_1','Fy_2','Fy_3','Fy_4'}, unit='N';
    case {'T_1','T_2','T_3','T_4','wheel_accel_residual_1','wheel_accel_residual_2','wheel_accel_residual_3','wheel_accel_residual_4'}, unit='N m';
    case {'omega_1','omega_2','omega_3','omega_4'}, unit='rad/s';
    case {'alpha_1','alpha_2','alpha_3','alpha_4','gamma_1','gamma_2','gamma_3','gamma_4'}, unit='deg';
    case {'kappa_1','kappa_2','kappa_3','kappa_4'}, unit='1';
    case {'wheel_accel_residual_1','wheel_accel_residual_2','wheel_accel_residual_3','wheel_accel_residual_4'}, unit='N m';
end
end

function channels=addDiagnosticChannels(channels,d,axisId,n)
spec={{'qss_status_code','QSS point status code','1',statusCodes(d),true,'none'}, ...
    {'qss_acceleration_residual_mps2','QSS acceleration residual','m/s^2',d.acceleration_residual_mps2,false,'linear'}, ...
    {'qss_yaw_acceleration_residual_radps2','QSS yaw acceleration residual','rad/s^2',d.yaw_acceleration_residual_radps2,false,'linear'}, ...
    {'qss_torque_residual_Nm','QSS wheel torque residual','N m',d.torque_residual_Nm,false,'linear'}, ...
    {'qss_min_virtual_load_N','QSS minimum virtual load','N',d.min_virtual_load_N,false,'linear'}, ...
    {'qss_aero_residual_in','QSS aero residual','in',d.aero_residual_in,false,'linear'}, ...
    {'qss_eval_count','QSS solver evaluation count','1',d.eval_count,false,'linear'}};
for i=1:numel(spec)
    s=spec{i}; ch=makeDetailChannel(s{1},s{2},s{3},axisId,n); ch.values=double(s{4}(:));
    if s{5}, ch.valid=true(n,1); else, ch.valid=isfinite(ch.values); ch.values(~ch.valid)=NaN; end
    ch.interpolation=s{6}; ch.description='Independent per-point reconstruction diagnostic'; channels(end+1)=ch; %#ok<AGROW>
end
end

function codes=statusCodes(d)
codes=zeros(numel(d.status),1); codes(strcmp(d.status,'solved'))=1; codes(strcmp(d.status,'failed'))=2;
end

function value=fieldOr(s,name,default)
if isstruct(s) && isfield(s,name), value=s.(name); else, value=default; end
end

function c=emptyChannel()
c=struct('id','','label','','unit','','axis_id','','origin','qss_reconstructed','interpolation','none', ...
    'description','','coordinate_frame','vehicle','sign_convention','','values',zeros(0,1),'valid',false(0,1));
end
