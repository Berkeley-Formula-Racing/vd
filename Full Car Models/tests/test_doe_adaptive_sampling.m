%% test_doe_adaptive_sampling
% A removed surrogate or acquisition branch must prevent deterministic,
% diverse adaptive selections from being returned.
testDir = fileparts(mfilename('fullpath'));
modelRoot = fileparts(testDir);
cd(modelRoot); setup_paths

U = [0.10 0.10; 0.15 0.80; 0.25 0.35; 0.35 0.65; ...
     0.45 0.20; 0.55 0.90; 0.65 0.45; 0.75 0.15; ...
     0.80 0.70; 0.90 0.30; 0.30 0.95; 0.95 0.95];
M = table(true(12,1),U(:,1)+U(:,2),U(:,1)-U(:,2), ...
    ones(12,1),500-(U(:,1)-0.8).^2*100, ...
    'VariableNames',{'valid','response_a','response_b','constant_response', ...
    'objective_score'});

models = doeFitSurrogates(U,M,{'response_a','response_b','constant_response'});
assert(isfield(models.responses,'response_a'))
assert(isfield(models.responses,'response_b'))
assert(ismember("constant_response",string(models.skippedResponses)))
assert(isnumeric(models.feasibility) && isscalar(models.feasibility) && ...
    models.feasibility == 1)

ei = doeExpectedImprovement([5;3],[0;2],4);
assert(abs(ei(1)-1) < 1e-12)
assert(ei(2) > 0)

study = DOEStudyConfig();
study.adaptive.responses = {'response_a','response_b'};
study.batchSize = 4;
study.adaptive.minimumDistance = 0.03;
state = struct('U',U,'metricTable',M);
C = [0.05 0.95;0.95 0.05;0.80 0.50;0.50 0.80; ...
     0.20 0.20;0.70 0.70;0.90 0.90];

study.mode = "sensitivity";
[Us,Ss] = doeSelectAdaptiveBatch(state,study,C);
assert(size(Us,1) == 4 && all(Ss.source == "sensitivity"))
assert(all(Ss.nearest_existing_distance >= study.adaptive.minimumDistance))

study.mode = "optimization";
[Uo,So] = doeSelectAdaptiveBatch(state,study,C);
assert(size(Uo,1) == 4 && all(So.source == "optimization"))

study.mode = "hybrid";
study.adaptive.hybridSensitivityFraction = 0.5;
[Uh,Sh] = doeSelectAdaptiveBatch(state,study,C);
assert(sum(Sh.source == "sensitivity") == 2)
assert(sum(Sh.source == "optimization") == 2)
assert(size(unique(Uh,'rows'),1) == 4)
assert(all(ismember(Sh.mode,"hybrid")))
