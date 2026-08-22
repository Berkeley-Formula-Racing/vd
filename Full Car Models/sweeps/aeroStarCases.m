function [carCell,plan] = aeroStarCases(baseCell,levels)
% One-at-a-time ("star") case set for sweeping ClA, CdA and CoP together.
%
% carConfig builds the CROSS PRODUCT of every vector it is given, so making
% all three aero parameters 5-level vectors there costs 125 full g-g solves.
% Nothing downstream uses the interactions: aeroSensitivity only ever
% differences cases that match the baseline in the other two parameters, so
% 112 of those 125 are computed and thrown away. This builds the 13 that are
% not -- one baseline plus four perturbations along each axis.
%
% inputs  baseCell: carConfig() output. Must be a SINGLE case (leave every
%                   aeroParams entry scalar in carConfig -- this function
%                   does the sweeping now).
%         levels:   optional struct
%                     .ClA  multipliers on the baseline ClA
%                           (default [0.90 0.95 1.05 1.10])
%                     .CdA  multipliers on the baseline CdA (same default)
%                     .CoP  ADDITIVE offsets in front downforce fraction,
%                           because CoP is already a fraction
%                           (default [-0.04 -0.02 +0.02 +0.04])
%                   A field set to [] skips that parameter entirely.
% outputs carCell:  n-by-2 cell in carConfig's own layout, case 1 = baseline
%         plan:     table describing each case (param, level, ClA/CdA/CoP)
%
% The baseline lands at the MIDDLE of every parameter's set of values, which
% is what aeroSensitivity's pickBaseline looks for, so it selects case 1 on
% its own with no opts.baseIdx needed.
%
% The accel-configuration car (column 2) is deliberately left at baseline.
% carConfig specifies the accel package separately (accel_cda / accel_cla) --
% it is a different wing setting, not a scaled version of this one. Accel
% times therefore come out identical across every case, which is expected and
% not a sign the sweep failed to apply.

if nargin < 2, levels = struct(); end
if size(baseCell,1) ~= 1
    error('aeroStarCases:notSingle', ...
        ['baseCell must hold exactly one case but holds %d. Set every ' ...
         'aeroParams entry in carConfig back to a scalar -- this function ' ...
         'builds the sweep, so a vector there multiplies the case count ' ...
         'again on top of it.'], size(baseCell,1));
end

mulDefault = [0.90 0.95 1.05 1.10];
lv.ClA = getOr(levels,'ClA',mulDefault);
lv.CdA = getOr(levels,'CdA',mulDefault);
lv.CoP = getOr(levels,'CoP',[-0.04 -0.02 0.02 0.04]);

base     = baseCell{1,1};
baseAcc  = baseCell{1,2};
a0       = base.aero;

% baseline must not appear twice: a multiplier of exactly 1 (or a CoP offset
% of 0) would duplicate case 1 and hand aeroSensitivity a zero dp to divide by
lv.ClA = lv.ClA(abs(lv.ClA-1) > 1e-12);
lv.CdA = lv.CdA(abs(lv.CdA-1) > 1e-12);
lv.CoP = lv.CoP(abs(lv.CoP)   > 1e-12);

vals = {a0.cla*lv.ClA(:).', a0.cda*lv.CdA(:).', a0.D_f + lv.CoP(:).'};
names = {'ClA','CdA','CoP'};

if any(vals{3} <= 0 | vals{3} >= 1)
    error('aeroStarCases:badCoP', ...
        ['a CoP offset puts the front downforce fraction outside (0,1): ' ...
         'baseline %.3f with offsets [%s]'], a0.D_f, num2str(lv.CoP));
end

% case 1 is the baseline, then each parameter's perturbations in turn
cla = a0.cla; cda = a0.cda; cop = a0.D_f;
param = {'baseline'}; level = 0;
for p = 1:3
    for k = 1:numel(vals{p})
        cla(end+1) = a0.cla; cda(end+1) = a0.cda; cop(end+1) = a0.D_f; %#ok<AGROW>
        switch p
            case 1, cla(end) = vals{p}(k);
            case 2, cda(end) = vals{p}(k);
            case 3, cop(end) = vals{p}(k);
        end
        param{end+1} = names{p};       %#ok<AGROW>
        level(end+1) = vals{p}(k);     %#ok<AGROW>
    end
end

n = numel(cla);
carCell = cell(n,2);
for i = 1:n
    c = base;
    c.aero = Aero(cda(i),cla(i),cop(i),a0.cla_p_deg_p,a0.D_p_deg_p);
    carCell{i,1} = c;
    carCell{i,2} = baseAcc;
end

plan = table((1:n)',string(param(:)),level(:),cla(:),cda(:),cop(:), ...
    'VariableNames',{'case','param','level','ClA','CdA','CoP'});
end

function v = getOr(s,f,d)
if isfield(s,f), v = s.(f); else, v = d; end
end
