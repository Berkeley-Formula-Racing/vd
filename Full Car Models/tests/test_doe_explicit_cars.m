%% test_doe_explicit_cars
testDir = fileparts(mfilename('fullpath'));
modelRoot = fileparts(testDir);
cd(modelRoot); setup_paths

[base,eventParams,X0,B] = carConfig();
assert(size(base,1)==1 && size(base,2)==2)
assert(height(X0)==1 && height(B)==1)
assert(abs(base{1,1}.M-226) < 1e-10)
assert(abs(base{1,1}.W_b-62*0.0254) < 1e-12)
assert(abs(base{1,1}.aero.cla-3.969) < 1e-12)
assert(isfield(eventParams,'winning_time'))

Q = table([155;170],[0.30;0.50],[3.5;4.2], ...
    'VariableNames',{'mass','R_sf','cla'});
[cars,~,X] = carConfig("Explicit",Q);
assert(size(cars,1)==2)
assert(abs(cars{1,1}.M-(155+64)) < 1e-10)
assert(abs(cars{2,1}.R_sf-0.50) < 1e-12)
assert(abs(cars{1,1}.aero.cla-3.5) < 1e-12)
assert(all(abs(X.p_i-12) < 1e-12))

[lhsCars,~,lhsX] = carConfig("LHS",3);
[~,warningId] = lastwarn;
assert(strcmp(warningId,'carConfig:legacyDOE'))
assert(size(lhsCars,1)==3 && height(lhsX)==3)
assert(all(lhsX.mass >= 162*0.95 & lhsX.mass <= 162*1.05))
assert(any(abs(lhsX.mass-162) < 1e-12))

study = DOEStudyConfig();
resolved = doeResolveStudy(study,B);
rng(19); expectedNext = rand;
rng(19); [U,T] = doeInitialDesign(resolved,3,7);
assert(abs(rand-expectedNext) < 1e-12)
[U2,T2] = doeInitialDesign(resolved,3,7);
assert(isequal(U,U2) && isequal(T,T2))
baselineU = (resolved.parameters.baseline-resolved.parameters.lowerPhysical) ./ ...
    (resolved.parameters.upperPhysical-resolved.parameters.lowerPhysical);
assert(any(all(abs(U-baselineU') < 1e-12,2)))
assert(size(U,1)==3 && height(T)==3)
[randomU,randomT] = doeInitialDesign(resolved,3,7,"Random");
assert(size(randomU,1)==3 && height(randomT)==3)
