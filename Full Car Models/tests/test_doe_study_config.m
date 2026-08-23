%% test_doe_study_config
addpath(fileparts(fileparts(mfilename('fullpath'))))
setup_paths
baseline = table(162,0.34,-1,0,'VariableNames', ...
    {'mass','R_sf','gamma_f','static_r_toe'});
study = struct();
study.parameters = table( ...
    ["mass";"R_sf";"gamma_f";"static_r_toe"], ...
    [-5;-10;-2;-0.2],[5;10;0;0.2], ...
    ["percent";"percent";"absolute";"absolute"], ...
    'VariableNames',{'name','lower','upper','rangeType'});
R = doeResolveStudy(study,baseline);
assert(abs(R.parameters.lowerPhysical(1)-153.9) < 1e-12)
assert(abs(R.parameters.upperPhysical(2)-0.374) < 1e-12)
T = R.toPhysical([0 0 0 0; 1 1 1 1]);
assert(isequal(T.Properties.VariableNames,cellstr(study.parameters.name)'))
assert(abs(T.gamma_f(1)+2) < 1e-12 && abs(T.gamma_f(2)) < 1e-12)
assert(R.signature == doeResolveStudy(study,baseline).signature)

bad = study; bad.parameters.rangeType(3) = "percent";
assertError(@() doeResolveStudy(bad,baseline),'doeResolveStudy:badPercentBaseline')
bad = study; bad.parameters.name(1) = "not_a_car_parameter";
assertError(@() doeResolveStudy(bad,baseline),'doeResolveStudy:unknownParameter')

function assertError(f,id)
try, f(); error('test:missingError','Expected %s',id)
catch ME, assert(strcmp(ME.identifier,id),ME.message)
end
end
