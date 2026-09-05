function [data,manifest] = lc0NDLoadFreeRolling(cfg)
%LC0NDLOADFREEROLLING Load configured 16x7.5x10 LC0 TTC runs in SI units.
%   TTC stores FZ as negative under compression and force in lbf. This
%   boundary exposes positive normal load and forces in N without changing
%   the measured FY sign convention.

if nargin < 1 || isempty(cfg), cfg = lc0NDConfig(); end
validateConfig(cfg);
target = cfg.target;
required = {'SA','SL','FY','FZ','IA','P','V'};
lbfToN = 4.4482216152605;
mphToMps = 0.44704;
fields = {'slipAngle_deg','slipRatio','Fy_N','Fz_N','camber_deg', ...
    'pressure_psi','speed_mps','runId'};
data = cell2struct(cell(size(fields)),fields,2);
manifest = table(zeros(0,1),strings(0,1),zeros(0,1),zeros(0,1), ...
    'VariableNames',{'run_id','file','n_samples','max_abs_slip'});

for runId = target.runIds(:)'
    file = fullfile(target.runDirectory,sprintf('%s%d.mat',target.filePrefix,runId));
    if ~isfile(file)
        error('lc0NDLoadFreeRolling:missingRun','Configured TTC run is missing: %s',file)
    end
    S = load(file);
    missing = setdiff(required,fieldnames(S));
    if ~isempty(missing)
        error('lc0NDLoadFreeRolling:missingChannel', ...
            'Run %d lacks TTC channel(s): %s',runId,strjoin(missing,', '))
    end
    n = numel(S.SA);
    for name = required
        value = S.(name{1});
        if ~isnumeric(value) || numel(value) ~= n || any(~isfinite(value(:)))
            error('lc0NDLoadFreeRolling:badChannel', ...
                'Run %d channel %s must be a finite vector with %d samples.', ...
                runId,name{1},n)
        end
    end
    slip = S.SL(:);
    maxSlip = max(abs(slip));
    if maxSlip > target.maxAbsSlip
        error('lc0NDLoadFreeRolling:notFreeRolling', ...
            'Run %d has max |SL|=%g, exceeding configured free-rolling limit %g.', ...
            runId,maxSlip,target.maxAbsSlip)
    end

    data.slipAngle_deg = [data.slipAngle_deg; S.SA(:)]; %#ok<AGROW>
    data.slipRatio = [data.slipRatio; slip]; %#ok<AGROW>
    data.Fy_N = [data.Fy_N; S.FY(:)*lbfToN]; %#ok<AGROW>
    data.Fz_N = [data.Fz_N; -S.FZ(:)*lbfToN]; %#ok<AGROW>
    data.camber_deg = [data.camber_deg; S.IA(:)]; %#ok<AGROW>
    data.pressure_psi = [data.pressure_psi; S.P(:)]; %#ok<AGROW>
    data.speed_mps = [data.speed_mps; S.V(:)*mphToMps]; %#ok<AGROW>
    data.runId = [data.runId; repmat(runId,n,1)]; %#ok<AGROW>
    manifest = [manifest; table(runId,string(file),n,maxSlip, ...
        'VariableNames',manifest.Properties.VariableNames)]; %#ok<AGROW>
end

if any(data.Fz_N <= 0)
    error('lc0NDLoadFreeRolling:nonpositiveLoad', ...
        'Configured free-rolling runs contain nonpositive normal loads.')
end
end

function validateConfig(cfg)
if ~isstruct(cfg) || ~isscalar(cfg) || ~isfield(cfg,'target') || ...
        ~isstruct(cfg.target) || ~isscalar(cfg.target)
    error('lc0NDLoadFreeRolling:badConfig','cfg.target must be a scalar struct.')
end
required = {'runDirectory','runIds','filePrefix','maxAbsSlip'};
if ~all(isfield(cfg.target,required))
    error('lc0NDLoadFreeRolling:badConfig', ...
        'cfg.target requires runDirectory, runIds, filePrefix, and maxAbsSlip.')
end
if ~isfolder(cfg.target.runDirectory) || ~isnumeric(cfg.target.runIds) || ...
        isempty(cfg.target.runIds) || any(cfg.target.runIds ~= floor(cfg.target.runIds)) || ...
        ~isnumeric(cfg.target.maxAbsSlip) || ~isscalar(cfg.target.maxAbsSlip) || ...
        ~isfinite(cfg.target.maxAbsSlip) || cfg.target.maxAbsSlip < 0
    error('lc0NDLoadFreeRolling:badConfig','Invalid free-rolling TTC configuration.')
end
end
