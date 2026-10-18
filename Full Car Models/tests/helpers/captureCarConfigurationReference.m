function referencePath = captureCarConfigurationReference(artifactDir)
%CAPTURECARCONFIGURATIONREFERENCE Save before-change construction evidence.
if nargin < 1 || isempty(artifactDir)
    helperDir = fileparts(mfilename('fullpath'));
    modelRoot = fileparts(fileparts(helperDir));
    repoRoot = fileparts(modelRoot);
    artifactDir = fullfile(repoRoot,'artifacts','car_configuration_refactor', ...
        '2026-10-04-reference');
end
if ~isfolder(artifactDir), mkdir(artifactDir); end
referencePath = fullfile(artifactDir,'construction-reference.mat');
if isfile(referencePath)
    error('captureCarConfigurationReference:exists', ...
        'Refusing to overwrite existing reference %s.',referencePath);
end

helperDir = fileparts(mfilename('fullpath'));
modelRoot = fileparts(fileparts(helperDir));
repoRoot = fileparts(modelRoot);
run(fullfile(modelRoot,'entrypoints','bootstrap.m'));
addpath(helperDir);

reference = struct();
reference.metadata = captureMetadata(repoRoot,modelRoot);
reference.metadata.artifactDirectory = artifactDir;
reference.metadata.capturePhase = 'construction-before-registry';
fprintf('Baseline construction capture in %s\n',artifactDir);

[cars,eventParams,designTable,baselineTable] = carConfig();
reference.qss.default.carCell = cars;
reference.qss.default.lapCar = carConfigurationSnapshot(cars{1,1});
reference.qss.default.accelCar = carConfigurationSnapshot(cars{1,2});
reference.qss.default.designTable = designTable;
reference.qss.default.designTableVariableNames = designTable.Properties.VariableNames;
reference.qss.default.baselineTable = baselineTable;
reference.qss.default.eventParams = eventParams;
reference.qss.default.torqueProbes = torqueProbes(cars{1,1});
fprintf('Captured scalar QSS roles and %d ordered design columns.\n',width(designTable));

explicitInput = table([155;170],[0.30;0.50],[3.5;4.2], ...
    'VariableNames',{'mass','R_sf','cla'});
[explicitCars,explicitEvents,explicitDesign,explicitBaseline] = ...
    carConfig("Explicit",explicitInput);
reference.qss.explicit.input = explicitInput;
reference.qss.explicit.carCell = explicitCars;
reference.qss.explicit.lapCars = {carConfigurationSnapshot(explicitCars{1,1}); ...
    carConfigurationSnapshot(explicitCars{2,1})};
reference.qss.explicit.accelCars = {carConfigurationSnapshot(explicitCars{1,2}); ...
    carConfigurationSnapshot(explicitCars{2,2})};
reference.qss.explicit.eventParams = explicitEvents;
reference.qss.explicit.designTable = explicitDesign;
reference.qss.explicit.baselineTable = explicitBaseline;
fprintf('Captured two-row Explicit QSS design.\n');

[ rampBaselineCar,rampConfig ] = carConfigBaseline();
rampAssets = rampSpeed.loadRampAssets(rampConfig,rampConfig.defaultSetup);
[rampBuiltCar,rampDerived] = buildSingleRampCar( ...
    rampConfig.fixedParameters,rampConfig.defaultSetup,rampAssets.map,rampAssets);
reference.ramp.config = rampConfig;
reference.ramp.carConfigBaselineCar = rampBaselineCar;
reference.ramp.builderCar = rampBuiltCar;
reference.ramp.snapshot = carConfigurationSnapshot(rampBuiltCar);
reference.ramp.derived = rampDerived;
reference.ramp.assetFingerprints = rampAssets.fingerprints;
reference.ramp.defaultMapPath = rampConfig.defaultMapPath;
reference.ramp.torqueProbes = torqueProbes(rampBuiltCar);
fprintf('Captured frozen ramp setup, car, derived fields, and asset identity.\n');

models = ["legacy","lc0_nd"];
uncertainties = ["low","nominal","high"];
selectionInput = table(162,'VariableNames',{'mass'});
for modelIndex = 1:numel(models)
    for uncertaintyIndex = 1:numel(uncertainties)
        model = models(modelIndex);
        uncertainty = uncertainties(uncertaintyIndex);
        [selectionCars,~,selectionDesign] = ...
            carConfig("Explicit",selectionInput,model,uncertainty);
        key = char(model + "_" + uncertainty);
        item = struct();
        item.model = model;
        item.uncertainty = uncertainty;
        item.designTable = selectionDesign;
        item.carCell = selectionCars;
        item.lapCar = carConfigurationSnapshot(selectionCars{1,1});
        item.accelCar = carConfigurationSnapshot(selectionCars{1,2});
        item.tireEvaluation = tireEvaluation(selectionCars{1,1}.tire);
        reference.qss.tireSelections.(key) = item;
        fprintf('Captured %s/%s construction and tire probe.\n',model,uncertainty);
    end
end

save(referencePath,'reference','-v7.3');
fprintf('Saved before-change construction evidence: %s\n',referencePath);
end

function metadata = captureMetadata(repoRoot,modelRoot)
metadata = struct();
metadata.repositoryRoot = repoRoot;
metadata.modelRoot = modelRoot;
metadata.workingDirectory = pwd;
metadata.matlabVersion = version;
metadata.matlabRelease = version('-release');
metadata.creationTime = datetime('now');
metadata.gitCommit = strtrim(systemText( ...
    'git -c safe.directory=C:/Users/johny/.codex/worktrees/car-configuration-refactor/VD rev-parse HEAD',repoRoot));
metadata.gitBranch = strtrim(systemText( ...
    'git -c safe.directory=C:/Users/johny/.codex/worktrees/car-configuration-refactor/VD branch --show-current',repoRoot));
metadata.gitStatus = systemText( ...
    'git -c safe.directory=C:/Users/johny/.codex/worktrees/car-configuration-refactor/VD status --short --branch',repoRoot);
functionNames = {'Car','Aero','Powertrain','Tire2','NondimensionalTire','AeroMap', ...
    'Events2','carConfig','carConfigBaseline','parameters_loop', ...
    'buildCarFromParameterSet','buildSingleRampCar','doeResolveStudy', ...
    'gg2','makeGG','rampSpeed.runLateralRamp','rampSpeed.runLongitudinalRamp', ...
    'buildTelemetryResult','exportQSSResults'};
metadata.resolvedPaths = struct();
for k = 1:numel(functionNames)
    name = functionNames{k};
    key = matlab.lang.makeValidName(strrep(name,'.','_'));
    metadata.resolvedPaths.(key) = which(name);
end
installed = ver;
metadata.installedToolboxes = {installed.Name};
metadata.optimizationToolbox = license('test','Optimization_Toolbox');
metadata.parallelComputingToolbox = license('test','Distrib_Computing_Toolbox');
metadata.statisticsToolbox = license('test','Statistics_Toolbox');
end

function output = systemText(command,workingDirectory)
old = pwd;
cleaner = onCleanup(@() cd(old)); %#ok<NASGU>
cd(workingDirectory);
[status,output] = system(command);
if status ~= 0
    output = sprintf('exit=%d\n%s',status,output);
end
end

function probes = torqueProbes(car)
curve = car.powertrain.torque_fn;
probes = struct('propertyClass',class(curve),'propertySize',size(curve), ...
    'rpm',[0 2500 5000 7500 10000 11500],'torque',[]);
if isa(curve,'function_handle')
    probes.functionText = func2str(curve);
    try
        probes.torque = arrayfun(curve,probes.rpm);
    catch ME
        probes.error = struct('identifier',ME.identifier,'message',ME.message);
    end
elseif isnumeric(curve) && size(curve,1) == 2
    probes.torque = interp1(curve(1,:),curve(2,:),probes.rpm,'linear','extrap');
else
    probes.value = curve;
end
end

function result = tireEvaluation(tire)
alpha = [0;1];
kappa = [0;0.02];
loadN = [500;500];
gammaDeg = [0;-1];
result.inputs = struct('alpha_deg',alpha,'kappa',kappa, ...
    'normalLoad_N',loadN,'camber_deg',gammaDeg);
result.Fx = tire.F_x(alpha,kappa,loadN,gammaDeg);
result.Fy = tire.F_y(alpha,kappa,loadN,gammaDeg);
end
