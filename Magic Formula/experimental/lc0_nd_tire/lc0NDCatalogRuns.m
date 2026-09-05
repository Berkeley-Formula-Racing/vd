function catalog = lc0NDCatalogRuns(rootDirectory)
%LC0NDCATALOGRUNS Inventory TTC files by their recorded slip coverage.
%   The catalog is intended to screen donor candidates before fitting. It
%   classifies data from channels, never from a directory name.

if nargin < 1 || ~isfolder(rootDirectory)
    error('lc0NDCatalogRuns:badDirectory', ...
        'rootDirectory must name an existing TTC data directory.');
end
files = dir(fullfile(rootDirectory,'**','*.mat'));
catalog = table(strings(0,1),strings(0,1),zeros(0,1),zeros(0,1), ...
    zeros(0,1),zeros(0,1),zeros(0,1),false(0,1),false(0,1), ...
    'VariableNames',{'file','folder','n_samples','max_abs_slip_ratio', ...
    'max_abs_slip_angle_deg','normal_load_min_N','normal_load_max_N', ...
    'is_free_rolling','has_combined_coverage'});
lbfToN = 4.4482216152605;
for i = 1:numel(files)
    file = fullfile(files(i).folder,files(i).name);
    S = load(file);
    required = {'SA','SL','FZ'};
    if ~all(isfield(S,required)) || ~isnumeric(S.SA) || ...
            ~isnumeric(S.SL) || ~isnumeric(S.FZ)
        continue
    end
    n = numel(S.SA);
    if n == 0 || numel(S.SL) ~= n || numel(S.FZ) ~= n
        continue
    end
    sa = S.SA(:);
    sl = S.SL(:);
    fz = -S.FZ(:)*lbfToN;
    valid = isfinite(sa) & isfinite(sl) & isfinite(fz) & fz > 0;
    if ~any(valid)
        continue
    end
    maxSlip = max(abs(sl(valid)));
    maxSA = max(abs(sa(valid)));
    isFree = maxSlip <= 1e-3;
    hasCombined = maxSlip >= 0.02 && maxSA >= 1;
    catalog = [catalog; table(string(files(i).name),string(files(i).folder), ...
        sum(valid),maxSlip,maxSA,min(fz(valid)),max(fz(valid)),isFree,hasCombined, ...
        'VariableNames',catalog.Properties.VariableNames)]; %#ok<AGROW>
end
catalog = sortrows(catalog,{'has_combined_coverage','max_abs_slip_ratio'}, ...
    {'descend','descend'});
end
