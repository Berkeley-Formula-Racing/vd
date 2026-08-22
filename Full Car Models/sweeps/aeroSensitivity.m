function S = aeroSensitivity(carCell,opts)
% Discrete-resweep sensitivity of the QSS envelope (and lap time) to the
% three aero parameters: ClA, CdA and CoP (front downforce fraction).
%
% Every case is a full independent re-solve of the g-g -- no linearisation --
% so the result captures re-optimisation of the operating point, not just a
% local gradient.
%
% inputs  carCell: cell array from carConfig() whose cars have ALREADY been
%                  through gg2 + makeGG. Sweep aero by making
%                  aeroParams.cla / .cda / .distribution vectors in carConfig.
%         opts:    optional struct
%                    .baseIdx   case treated as baseline (default: middle
%                               value of each swept parameter, else case 1)
%                    .times     table of lap times per case, one row per
%                               case, from aeroLapTimes (optional)
%                    .thetaGrid, .gLatGrid, .gLongGrid  override grids
%
% output  S: struct holding
%            .cases     table of ClA/CdA/CoP per case + which param varies
%            .grids     the common axes: vCar, theta, gLat, gLong
%            .env       per-case envelopes on those grids (see aeroEnvelope)
%            .metrics   stacked ggMetrics table across all cases
%            .sens      sensitivities, S.sens.ClA.dr_dp etc.
%            .times     lap times per case, if supplied
%
% Sensitivities are reported two ways:
%   dr_dp      absolute, metric units per unit of the parameter
%   dr_dp_pct  normalised, % change in metric per 1% change in the parameter
%              (for CoP, which is already a fraction, the normalised form is
%              per 0.01 of balance -- see S.sens.CoP.note)

if nargin < 2, opts = struct(); end

% struct('traces',someCell) silently builds a 1-by-N STRUCT ARRAY, one element
% per cell, rather than a scalar struct holding a cell. Everything downstream
% then fails with an unrelated "Too many input arguments".
if numel(opts) ~= 1
    error('aeroSensitivity:optsNotScalar', ...
        ['opts must be a scalar struct but is %d-by-%d. If you built it with ' ...
         'struct(), wrap cell-valued fields in braces:\n' ...
         '    struct(''times'',times,''traces'',{traces})\n' ...
         'or assign field by field:\n' ...
         '    opts = struct(); opts.times = times; opts.traces = traces;'], ...
         size(opts,1),size(opts,2));
end

cars = carCell(:,1);
n = numel(cars);
if n < 2
    error('aeroSensitivity:tooFewCases', ...
        ['need at least 2 cases; sweep aeroParams.cla, .cda or .distribution ' ...
         'as a vector in carConfig.m']);
end

%% case table
ClA = zeros(n,1); CdA = zeros(n,1); CoP = zeros(n,1);
for i = 1:n
    ClA(i) = cars{i}.aero.cla;
    CdA(i) = cars{i}.aero.cda;
    CoP(i) = cars{i}.aero.D_f;
end
S.cases = table((1:n)',ClA,CdA,CoP,'VariableNames',{'case','ClA','CdA','CoP'});

% Guard: every case must differ ONLY in aero. carConfig builds the cross
% product of all swept parameters, so leaving e.g. wheelbase as a vector
% would have us differencing cars with different geometry and booking the
% result as an aero sensitivity.
checkNonAeroConstant(cars);

params = {'ClA','CdA','CoP'};
vals   = {ClA,CdA,CoP};
varies = cellfun(@(v) numel(unique(round(v,12)))>1, vals);
if ~any(varies)
    error('aeroSensitivity:noVariation', ...
        'no aero parameter varies across the supplied cars');
end

%% baseline
if isfield(opts,'baseIdx')
    baseIdx = opts.baseIdx;
else
    baseIdx = pickBaseline(vals,varies,n);
end
S.baseIdx = baseIdx;
S.cases.isBaseline = (1:n)' == baseIdx;

%% common grids
vAll = [];
for i = 1:n
    if ~isempty(cars{i}.ss_info)
        vAll = [vAll; unique(cars{i}.ss_info(:,6))]; %#ok<AGROW>
    end
end
grids.vCar  = unique(round(vAll,6));
grids.theta = getOr(opts,'thetaGrid',linspace(0,180,73));    % 2.5 deg steps
grids.gLat  = getOr(opts,'gLatGrid', linspace(0,2.0,41));
grids.gLong = getOr(opts,'gLongGrid',linspace(-2.5,1.5,41));
S.grids = grids;

%% envelopes and metric tables
S.env = cell(n,1);
metricTbls = cell(n,1);
for i = 1:n
    S.env{i} = aeroEnvelope(cars{i},grids);
    label = sprintf('case%02d ClA=%.3f CdA=%.3f CoP=%.3f',i,ClA(i),CdA(i),CoP(i));
    metricTbls{i} = ggMetrics(cars{i},label);
    metricTbls{i}.case = repmat(i,height(metricTbls{i}),1);
end
S.metrics = vertcat(metricTbls{:});

%% sensitivities, one parameter at a time
envFields = {'r_theta','accel_at_gLat','brake_at_gLat','gLat_at_gLong', ...
             'maxGLat','maxGLong','maxGBrake'};
base = S.env{baseIdx};

for p = find(varies)
    name = params{p};
    v    = vals{p};
    % cases that differ from baseline in THIS parameter only
    others = setdiff(1:3,p);
    sameOthers = true(n,1);
    for q = others
        sameOthers = sameOthers & abs(vals{q}-vals{q}(baseIdx)) < 1e-12;
    end
    sel = find(sameOthers & abs(v-v(baseIdx)) > 1e-12);
    if isempty(sel)
        warning('aeroSensitivity:noCleanCases', ...
            ['%s varies but never in isolation from the other aero ' ...
             'parameters -- skipping. Sweep one parameter at a time, or ' ...
             'set opts.baseIdx to a case that shares the others.'],name);
        continue
    end

    Sp = struct();
    Sp.cases   = sel;
    Sp.pBase   = v(baseIdx);
    Sp.pValues = v(sel);
    Sp.dp      = v(sel) - v(baseIdx);

    for f = 1:numel(envFields)
        fn = envFields{f};
        B  = base.(fn);
        D  = nan([size(B) numel(sel)]);
        for k = 1:numel(sel)
            D(:,:,k) = (S.env{sel(k)}.(fn) - B)./Sp.dp(k);
        end
        Sp.([fn '_dp'])       = mean(D,3,'omitnan');
        Sp.([fn '_dp_spread'])= max(D,[],3,'omitnan') - min(D,[],3,'omitnan');
        Sp.([fn '_dp_all'])   = D;
        % Reliability: if the perturbation levels disagree on the SIGN of the
        % slope, the underlying solves found different local optima and the
        % averaged gradient is meaningless. Near max velocity the accel
        % envelope collapses and this happens routinely.
        if size(D,3) > 1
            lo = min(D,[],3,'omitnan'); hi = max(D,[],3,'omitnan');
            Sp.([fn '_reliable']) = ~(lo < 0 & hi > 0) & isfinite(Sp.([fn '_dp']));
        else
            Sp.([fn '_reliable']) = isfinite(Sp.([fn '_dp']));
        end
        % normalised: % metric change per % parameter change
        if strcmp(name,'CoP')
            Sp.([fn '_dp_norm']) = Sp.([fn '_dp'])*0.01;   % per 0.01 balance
        else
            Sp.([fn '_dp_norm']) = Sp.([fn '_dp']).*(Sp.pBase./B)*0.01;
        end
    end

    if strcmp(name,'CoP')
        Sp.note = 'normalised columns are per 0.01 (1 percentage point) of front balance';
    else
        Sp.note = 'normalised columns are % metric change per 1% parameter change';
    end
    S.sens.(name) = Sp;
end

%% lap traces (for the time-at-condition kernels)
if isfield(opts,'traces') && ~isempty(opts.traces)
    S.traces = opts.traces;
end

%% lap times
if isfield(opts,'times') && ~isempty(opts.times)
    S.times = opts.times;
    % which event decides "fastest"; autocross by default
    S.fastestBy = getOr(opts,'fastestBy','autocross');
    tCols = setdiff(S.times.Properties.VariableNames,{'case'});
    for p = find(varies)
        name = params{p};
        if ~isfield(S.sens,name), continue, end
        sel = S.sens.(name).cases;
        dp  = S.sens.(name).dp;
        for c = 1:numel(tCols)
            tb = S.times.(tCols{c})(baseIdx);
            dt = (S.times.(tCols{c})(sel) - tb)./dp;
            S.sens.(name).(['d_' tCols{c} '_dp']) = mean(dt,'omitnan');
            S.sens.(name).(['d_' tCols{c} '_dp_all']) = dt;
        end

        lap  = S.times.(S.fastestBy);
        cand = [baseIdx; sel(:)];
        [~,k] = min(lap(cand));
        S.sens.(name).fastestCase = cand(k);

        [~,k2] = min(lap(sel));
        S.sens.(name).compareCase = sel(k2);
        S.sens.(name).baselineIsFastest = (cand(k) == baseIdx);
        S.sens.(name).compareDp = vals{p}(sel(k2)) - vals{p}(baseIdx);
    end
end
end

function checkNonAeroConstant(cars)
% errors if any non-aero vehicle parameter varies across the cases
carProps = {'M','W_b','l_f','t_f','t_r','h_g','h_rf','h_rr','R_sf','I_zz','R', ...
            'static_gamma_f','static_gamma_r','camber_compliance_f', ...
            'camber_compliance_r','ackermann','static_r_toe'};
ptProps  = {'final_drive','brake_distribution','max_braking_torque','redline', ...
            'shift_point','drivetrain_efficiency','G_d2_driving'};
tireProps = {'p_i','friction_scaling_factor'};

offenders = {};
offenders = [offenders, scan(cars,carProps,@(c) c)];
offenders = [offenders, scan(cars,ptProps,@(c) c.powertrain)];
offenders = [offenders, scan(cars,tireProps,@(c) c.tire)];

if ~isempty(offenders)
    error('aeroSensitivity:nonAeroVariation', ...
        ['these non-aero parameters vary across the cases: %s\n' ...
         'carConfig takes the cross product of every swept parameter, so an ' ...
         'aero sweep must be the ONLY vector in carConfig.m. Set the others ' ...
         'back to single values and re-run.'], strjoin(offenders,', '));
end
end

function bad = scan(cars,props,getter)
bad = {};
for k = 1:numel(props)
    v = zeros(numel(cars),1);
    ok = true;
    for i = 1:numel(cars)
        s = getter(cars{i});
        if ~isprop(s,props{k}) && ~isfield(s,props{k}), ok = false; break, end
        x = s.(props{k});
        if ~isscalar(x) || ~isnumeric(x), ok = false; break, end
        v(i) = x;
    end
    if ok && numel(unique(round(v,10))) > 1
        bad{end+1} = props{k};
    end
end
end

function b = pickBaseline(vals,varies,n)
% prefer the case sitting at the middle value of every swept parameter
cand = true(n,1);
for p = find(varies)
    u = unique(round(vals{p},12));
    mid = u(ceil(numel(u)/2));
    cand = cand & abs(vals{p}-mid) < 1e-12;
end
b = find(cand,1);
if isempty(b), b = 1; end
end

function v = getOr(s,f,d)
if isfield(s,f), v = s.(f); else, v = d; end
end
