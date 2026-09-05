function S = aeroMapSensitivity(carCell,plan,opts)
%AEROMAPSENSITIVITY Sensitivity of solved g-g envelopes to aeromap inputs.
%   PLAN is created by aeroMapStarCases. Each perturbation is a full re-solve,
%   so the reported derivatives include the coupled ride-height/aero response.

if nargin < 3, opts = struct(); end
if numel(opts) ~= 1
    error('aeroMapSensitivity:optsNotScalar','opts must be a scalar struct.');
end
if size(carCell,1) ~= height(plan)
    error('aeroMapSensitivity:caseCount', ...
        'carCell and plan must contain the same number of cases.');
end

info = aeroMapSensitivityPlan(plan);
cars = carCell(:,1);
if any(cellfun(@(car) isempty(car.ss_info),cars))
    error('aeroMapSensitivity:unsolvedCars', ...
        'every case must have solved g-g data; run gg2 and makeGG first.');
end

S = struct();
S.cases = plan;
S.baseIdx = info.baseIdx;
S.cases.isBaseline = false(height(plan),1);
S.cases.isBaseline(S.baseIdx) = true;
S.paramInfo = info.parameters;
S.grids = commonGrids(cars,opts);

S.env = cell(numel(cars),1);
metricTbls = cell(numel(cars),1);
for i = 1:numel(cars)
    S.env{i} = aeroEnvelope(cars{i},S.grids);
    label = sprintf('case%02d %s=%+.4g',plan.case(i), ...
        char(plan.swept_parameter(i)),plan.swept_value(i));
    metricTbls{i} = ggMetrics(cars{i},label);
    metricTbls{i}.case = repmat(plan.case(i),height(metricTbls{i}),1);
end
S.metrics = vertcat(metricTbls{:});

fields = {'r_theta','accel_at_gLat','brake_at_gLat','gLat_at_gLong', ...
    'maxGLat','maxGLong','maxGBrake'};
S.sens = struct();
for p = 1:numel(info.parameters)
    spec = info.parameters(p);
    name = char(spec.name);
    values = plan.(char(spec.column));
    selected = find(string(plan.swept_parameter) == spec.name);
    verifyIsolation(plan,info,p,selected,S.baseIdx);
    if isempty(selected) || any(abs(values(selected)-values(S.baseIdx)) < 1e-12)
        error('aeroMapSensitivity:zeroPerturbation', ...
            '%s must contain nonzero perturbations from the baseline.',name);
    end

    P = struct();
    P.cases = selected;
    P.pBase = values(S.baseIdx);
    P.pValues = values(selected);
    P.dp = values(selected)-P.pBase;
    P.label = spec.label;
    P.unit = spec.unit;
    P.note = sprintf('absolute derivative per %s of %s',spec.unit,spec.label);

    for f = 1:numel(fields)
        field = fields{f};
        baseline = S.env{S.baseIdx}.(field);
        D = nan([size(baseline),numel(selected)]);
        for k = 1:numel(selected)
            D(:,:,k) = (S.env{selected(k)}.(field)-baseline)./P.dp(k);
        end
        P.([field '_dp']) = mean(D,3,'omitnan');
        P.([field '_dp_spread']) = max(D,[],3,'omitnan') - min(D,[],3,'omitnan');
        P.([field '_dp_all']) = D;
        if size(D,3) > 1
            lo = min(D,[],3,'omitnan');
            hi = max(D,[],3,'omitnan');
            P.([field '_reliable']) = ~(lo < 0 & hi > 0) & isfinite(P.([field '_dp']));
        else
            P.([field '_reliable']) = isfinite(P.([field '_dp']));
        end
    end
    S.sens.(name) = P;
end

if isfield(opts,'traces') && ~isempty(opts.traces), S.traces = opts.traces; end
if isfield(opts,'times') && ~isempty(opts.times)
    S = addLapTimes(S,plan,opts.times,getOr(opts,'fastestBy','autocross'));
end
end

function grids = commonGrids(cars,opts)
vAll = [];
for i = 1:numel(cars)
    vAll = [vAll; unique(cars{i}.ss_info(:,6))]; %#ok<AGROW>
end
grids.vCar = unique(round(vAll,6));
grids.theta = getOr(opts,'thetaGrid',linspace(0,180,73));
grids.gLat = getOr(opts,'gLatGrid',linspace(0,2.0,41));
grids.gLong = getOr(opts,'gLongGrid',linspace(-2.5,1.5,41));
end

function verifyIsolation(plan,info,p,selected,baseIdx)
for q = 1:numel(info.parameters)
    if q == p, continue, end
    column = char(info.parameters(q).column);
    if any(abs(plan.(column)(selected)-plan.(column)(baseIdx)) > 1e-12)
        error('aeroMapSensitivity:notStarDesign', ...
            'cases for %s also change %s.',info.parameters(p).name,info.parameters(q).name);
    end
end
end

function S = addLapTimes(S,plan,times,fastestBy)
S.times = times;
S.fastestBy = fastestBy;
timeColumns = setdiff(times.Properties.VariableNames,{'case'});
for p = 1:numel(S.paramInfo)
    spec = S.paramInfo(p);
    name = char(spec.name);
    P = S.sens.(name);
    for c = 1:numel(timeColumns)
        column = timeColumns{c};
        baseline = times.(column)(S.baseIdx);
        delta = (times.(column)(P.cases)-baseline)./P.dp;
        P.(['d_' column '_dp']) = mean(delta,'omitnan');
        P.(['d_' column '_dp_all']) = delta;
    end
    lap = times.(fastestBy);
    candidates = [S.baseIdx; P.cases(:)];
    [~,index] = min(lap(candidates));
    P.fastestCase = candidates(index);
    [~,index] = min(lap(P.cases));
    P.compareCase = P.cases(index);
    P.baselineIsFastest = P.fastestCase == S.baseIdx;
    P.compareDp = P.pValues(index)-P.pBase;
    S.sens.(name) = P;
end
end

function value = getOr(s,name,default)
if isfield(s,name), value = s.(name); else, value = default; end
end
