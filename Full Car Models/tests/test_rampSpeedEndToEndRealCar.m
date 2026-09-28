function tests = test_rampSpeedEndToEndRealCar
tests = functiontests(localfunctions);
end

function testBaselineRunsBothRampTypesAtAcceptanceSpeeds(testCase)
ensureRampSpeedPath();
[baselineCar,config] = carConfigBaseline();
[cars,cases] = rampSpeed.buildSetupCatalog(config,config.defaultSetup);
verifySize(testCase,cars,[1 1]);
verifyClass(testCase,cars{1},"Car");
verifyEqual(testCase,string(cases(1).id),"ramp-baseline");

speeds = [5 10 15 17.5 20 22.5 25];
rampTypes = ["longitudinal","lateral"];
allowedStatuses = ["planned","running","converged","near_feasible", ...
    "infeasible","solver_failed","cancelled","complete","completed","failed"];

for rampIndex = 1:numel(rampTypes)
    request = makeRealCarRequest(rampTypes(rampIndex),speeds);
    [study,events] = rampSpeed.runStudy(cars,cases,request,struct());
    verifyNotEqual(testCase,string(study.status),"failed");
    verifyEqual(testCase,numel(study.runs),1);
    verifyGreaterThanOrEqual(testCase,numel(events),1);
    run = study.runs(1);
    verifyEqual(testCase,run.perSpeed.speed_mps,speeds(:),"AbsTol",0);
    verifyEqual(testCase,height(run.perSpeed),numel(speeds));
    verifyTrue(testCase,all(ismember(string(run.perSpeed.status), ...
        allowedStatuses)));
    verifyTrue(testCase,all(isfinite(run.perSpeed.speed_mps)));
    invalid = ~logical(run.perSpeed.valid);
    if any(invalid)
        if rampTypes(rampIndex) == "longitudinal"
            verifyTrue(testCase,all(isnan(run.perSpeed.aLong_mps2(invalid))));
        else
            verifyTrue(testCase,all(isnan( ...
                run.perSpeed.aLat_free_mps2(invalid))));
            verifyTrue(testCase,all(isnan( ...
                run.perSpeed.mechanical_balance_front(invalid))));
        end
    end
    if rampTypes(rampIndex) == "longitudinal"
        verifyEqual(testCase,string(run.runMeta.rampModel.modelKind), ...
            "rampSpeedLite");
        verifyEqual(testCase,string(run.runMeta.rampModel.powertrainModel), ...
            "continuousEnvelope");
        for diagnostic = run.raw.diagnostics.'
            if isfield(diagnostic.diagnostics,"envelopeSource")
                verifyEqual(testCase,string(diagnostic.diagnostics.envelopeSource), ...
                    "cachedRampModel");
            end
        end
    end
end

clear baselineCar config cars cases
end

function testLegacyWrapperNamesCanonicalBaselinePath(testCase)
ensureRampSpeedPath();
wrapperPath = which("runRampSpeedStudy");
verifyNotEmpty(testCase,wrapperPath);
source = fileread(wrapperPath);
verifyNotEmpty(testCase,regexp(source,"carConfigBaseline","once"));
verifyNotEmpty(testCase,regexp(source,"rampSpeed\.runStudy","once"));
verifyEmpty(testCase,regexp(source,"\bcarConfig\s*\(","once"));
verifyEmpty(testCase,regexp(source,"\brampSweep\s*\(","once"));
verifyEmpty(testCase,regexp(source,"\bparfor\b","once"));
end

function request = makeRealCarRequest(rampType,speeds)
settings = struct("speeds",double(speeds(:).'), ...
    "nRamp",4,"nBisect",0,"mode","coast","verbose",false, ...
    "solverProfile","fastPreview","speedGrid",struct("mode","fixed"));
request = struct("rampType",string(rampType),"settings",settings, ...
    "parallelRequested",false,"numWorkers",0,"checkpointPath","", ...
    "appVersion","ramp-speed-real-e2e","runCaseFcn",@rampSpeed.runCase);
end

function ensureRampSpeedPath()
root = fileparts(fileparts(mfilename("fullpath")));
addpath(genpath(root));
end
