function donor = lc0NDLoadDonor(cfg)
%LC0NDLOADDONOR Load one provisional combined-slip TTC donor in SI units.
%   Donor force signs retain TTC convention. Fz is exposed positive under
%   compression so force coefficients use Fx_N/Fz_N and Fy_N/Fz_N.

if nargin < 1 || ~isstruct(cfg) || ~isfield(cfg,'file') || ~isfile(cfg.file)
    error('lc0NDLoadDonor:badConfig','CFG.file must name an existing TTC MAT file.');
end
S = load(cfg.file);
required = {'SA','SL','FX','FY','FZ','IA','P','V'};
missing = setdiff(required,fieldnames(S));
if ~isempty(missing)
    error('lc0NDLoadDonor:missingChannel', ...
        'Donor file lacks TTC channel(s): %s.',strjoin(missing,', '));
end
n = numel(S.SA);
for i = 1:numel(required)
    value = S.(required{i});
    if ~isnumeric(value) || numel(value) ~= n || any(~isfinite(value(:)))
        error('lc0NDLoadDonor:badChannel', ...
            'TTC channel %s must be a finite vector with %d samples.',required{i},n);
    end
end
lbfToN = 4.4482216152605;
donor = struct('file',string(cfg.file), ...
    'slipAngle_deg',S.SA(:), ...
    'slipRatio',S.SL(:), ...
    'Fx_N',S.FX(:)*lbfToN, ...
    'Fy_N',S.FY(:)*lbfToN, ...
    'Fz_N',-S.FZ(:)*lbfToN, ...
    'camber_deg',S.IA(:), ...
    'pressure_psi',S.P(:), ...
    'speed_mps',S.V(:)*0.44704);
if any(donor.Fz_N <= 0)
    error('lc0NDLoadDonor:nonpositiveLoad', ...
        'Donor file contains nonpositive normal load.');
end
end
