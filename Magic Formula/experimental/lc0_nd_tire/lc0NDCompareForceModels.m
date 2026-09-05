function comparison = lc0NDCompareForceModels(model,cases,legacyEvaluator)
%LC0NDCOMPAREFORCEMODELS Evaluate experimental and legacy tyres at same states.

required = {'alpha_deg','slip_ratio','Fz_N','camber_deg'};
if ~istable(cases) || ~all(ismember(required,cases.Properties.VariableNames))
    error('lc0NDCompareForceModels:badCases', ...
        'CASES must be a table with alpha_deg, slip_ratio, Fz_N, camber_deg.');
end
if ~isa(legacyEvaluator,'function_handle')
    error('lc0NDCompareForceModels:badLegacyEvaluator', ...
        'legacyEvaluator must be a function handle returning [Fx,Fy].');
end
[experimentalFx,experimentalFy,info] = lc0NDEvaluate(model,cases.alpha_deg, ...
    cases.slip_ratio,cases.Fz_N);
n = height(cases);
legacyFx = nan(n,1); legacyFy = nan(n,1);
for i = 1:n
    [legacyFx(i),legacyFy(i)] = legacyEvaluator(cases.alpha_deg(i), ...
        cases.slip_ratio(i),cases.Fz_N(i),cases.camber_deg(i));
end
comparison = cases;
comparison.experimental_Fx_N = experimentalFx;
comparison.experimental_Fy_N = experimentalFy;
comparison.legacy_Fx_N = legacyFx;
comparison.legacy_Fy_N = legacyFy;
comparison.experimental_supported = info.is_supported;
comparison.experimental_utilization = info.utilization;
end
