function tests = test_qssHardening
% Regression tests for QSS validity and hot-path hardening.
tests = functiontests(localfunctions);
end

function setupOnce(~)
testDir = fileparts(mfilename('fullpath'));
addpath(fileparts(testDir));
setup_paths
end

function testConstraintUsesConfiguredPowertrainRedline(testCase)
[cars,~] = carConfig();
car = cars{1,1};

[P,engineRpm] = stateBetweenConfiguredAndLegacyRedline(car);
[c,~] = car.constraint1(P);

verifyGreaterThan(testCase,engineRpm,car.powertrain.redline);
verifyLessThan(testCase,engineRpm,13000);
verifyGreaterThan(testCase,c(1),0);
end

function testParametersLoopPassesPreloadedCamberModelData(testCase)
[cars,~] = carConfig();

verifyTrue(testCase,isfield(cars{1,1}.camberKinematics,'camberModelData'));
verifyTrue(testCase,isstruct(cars{1,1}.camberKinematics.camberModelData));
end

function testEvaluatorCarriesAndReusesRideHeightContext(testCase)
[cars,~] = carConfig();
car = cars{1,1};
P = [4,0.15,20,0.2,0.45,0,0,0.02,0.02];

evaluator = steadyStateEvaluator(car);
first = evaluator.evaluate(P);
P2 = P;
P2(2) = P2(2) + 0.01;
second = evaluator.evaluate(P2);

verifyTrue(testCase,isfield(first,'rideHeightContext'));
verifyTrue(testCase,isfield(first,'ssInfo'));
verifyTrue(testCase,first.aeroConverged);
verifyFalse(testCase,first.aeroOutsideMap);
verifyTrue(testCase,second.rideHeightContext.usedWarmStart);
verifyEqual(testCase,evaluator.equationCalls(),2);
end

function testSteadyStateConstraintRejectsInvalidAeroState(testCase)
state = struct( ...
    'engineRpm',1000,'beta',0,'Fzvirtual',[1 1 1 1], ...
    'latAccel',0,'yawAccel',0,'wheelAccel',[0 0 0 0], ...
    'aeroConverged',false,'aeroOutsideMap',false);

[c,ceq] = steadyStateConstraint4(state,[0 0 20 0 0 0 0 0 0],0);

verifyTrue(testCase,any(c > 0));
verifyTrue(testCase,all(isfinite([c ceq])));
end

function testCachedConstraintRejectsMissingRedline(testCase)
state = struct( ...
    'engineRpm',1000,'beta',0,'Fzvirtual',[1 1 1 1], ...
    'latAccel',0,'yawAccel',0,'wheelAccel',[0 0 0 0], ...
    'aeroConverged',true,'aeroOutsideMap',false, ...
    'aeroCoverageValid',true,'aeroResidualIn',0, ...
    'aeroResidualToleranceIn',1e-9);

[c,~] = steadyStateConstraint1(state);

verifyGreaterThan(testCase,c(1),0);
end

function testLargeAeroResidualIsRejected(testCase)
state = struct( ...
    'engineRpm',1000,'beta',0,'Fzvirtual',[1 1 1 1], ...
    'latAccel',0,'yawAccel',0,'wheelAccel',[0 0 0 0], ...
    'aeroConverged',true,'aeroOutsideMap',false, ...
    'aeroCoverageValid',true,'aeroResidualIn',1e-3, ...
    'aeroResidualToleranceIn',1e-9);

[c,~] = steadyStateConstraint4(state,[0 0 20 0 0 0 0 0 0],0);

verifyGreaterThan(testCase,max(c),0);
end

function testNearZeroYawRadiusConstraintRemainsFinite(testCase)
[cars,~] = carConfig();
car = cars{1,1};
for yawRate = [0,1e-12,-1e-12]
    P = [0,0,20,0,yawRate,0,0,0,0];
    [c,ceq] = car.constraint3(P,8.5);
    verifyTrue(testCase,all(isfinite([c ceq])));
    verifyEqual(testCase,ceq(1),P(3)-8.5*P(5),'AbsTol',1e-12);
end
end

function testCandidateValidityRequiresFiniteResidualAndSuccessfulExit(testCase)
verifyTrue(testCase,qssCandidateValidity(1,1e-5,1e-4));
verifyTrue(testCase,qssCandidateValidity(2,1e-5,1e-4));
verifyFalse(testCase,qssCandidateValidity(0,1e-5,1e-4));
verifyFalse(testCase,qssCandidateValidity(1,1e-2,1e-4));
verifyFalse(testCase,qssCandidateValidity(1,Inf,1e-4));
verifyEqual(testCase,qssConstraintResidual([-10;-1],0),0);
verifyEqual(testCase,qssConstraintResidual([0.2;-1],0),0.2);
end

function testMakeGGPreservesInvalidCellMask(testCase)
p = ParamSet(struct(),10);
p.maxLatFlag = 1;
p.maxLatValid = false;
p.maxLatx = zeros(1,9);
p.maxLatx(6) = 10;
p.maxLongFlag = 1;
p.maxLongValid = true;
p.maxLongLongAccel = 1;
p.maxLongLatAccel = 2;
p.maxLongAccelx = zeros(1,9);
p.maxBrakeFlag = 1;
p.maxBrakeValid = true;
p.maxBrakeLongDecel = -1;
p.maxBrakeLatAccel = 2;
p.maxBrakeDecelx = zeros(1,9);

result = makeGG(p,struct());

verifyFalse(testCase,result.ggMask.lateral);
verifyFalse(testCase,result.ggMask.acceleration);
verifyFalse(testCase,result.ggMask.braking);
verifyEmpty(testCase,result.ggPoints);
end

function testMakeGGHandlesBlankMissingParamSet(testCase)
p = ParamSet(struct(),10);
result = makeGG(p,struct());

verifyFalse(testCase,result.ggMask.lateral);
verifyFalse(testCase,result.ggMask.acceleration);
verifyFalse(testCase,result.ggMask.braking);
end

function testEnvelopeLookupPreservesMissingVelocityGaps(testCase)
accel = [1.0,0.5,20; 1.2,0.6,20; 1.1,0.4,30; 1.3,0.7,30];
brake = [-1.0,0.5,20; -1.2,0.6,20; -1.1,0.4,30; -1.3,0.7,30];
ggMask = struct( ...
    'velocity',[10;10;20;20;30;30], ...
    'acceleration',[false;false;true;true;true;true], ...
    'braking',[false;false;true;true;true;true]);

[F_accel,F_braking] = create_scattered_interpolants2(accel,brake,ggMask);

verifyTrue(testCase,isnan(F_accel(0.5,10)));
verifyTrue(testCase,isnan(F_braking(0.5,15)));
verifyTrue(testCase,isfinite(F_accel(0.5,25)));
end

function testEnvelopeLookupRejectsMissingBoundaryRows(testCase)
accel = [1.0,0.5,20; 1.2,0.6,20; 1.1,0.4,30; 1.3,0.7,30];
brake = [-1.0,0.5,20; -1.2,0.6,20; -1.1,0.4,30; -1.3,0.7,30];
ggMask = struct( ...
    'velocity',[10;10;20;20;30;30], ...
    'acceleration',[true;true;true;true;false;false], ...
    'braking',[true;true;true;true;false;false]);

[F_accel,F_braking] = create_scattered_interpolants2(accel,brake,ggMask);

verifyTrue(testCase,isnan(F_accel(0.5,25)));
verifyTrue(testCase,isnan(F_accel(0.5,35)));
verifyTrue(testCase,isnan(F_braking(0.5,25)));
end

function testCornerVelocityLookupPreservesMissingGgRows(testCase)
car = struct();
car.max_vel = 30;
car.ss_info = [ ...
    0,0,10,0,0,10; 0,0,10,0,0,10; ...
    0,0,30,0,0,30; 0,0,30,0,0,30];
car.ggMask = struct('velocity',[10,10;20,20;30,30], ...
    'lateral',[true,true;false,false;true,true]);

[~,radius,maxVelocity] = vel_cornering_from_gg(car);
gapIndex = find(abs(radius-12) < 1e-12,1);

verifyTrue(testCase,isnan(maxVelocity(gapIndex)));
end

function testCornerVelocityLookupHonoursGgVelocityCap(testCase)
car = struct();
car.max_vel = 33;
car.ss_info = [0,0,15,0,0,5; 0,0,15,0,0,5; ...
    0,0,15,0,0,28; 0,0,15,0,0,28];
car.ggMask = struct('velocity',[5,5;28,28], ...
    'lateral',[true,true;true,true]);

[~,~,maxVelocity] = vel_cornering_from_gg(car);

verifyLessThanOrEqual(testCase,max(maxVelocity),28 + 1e-12);
end

function testConstraint9IncludesLateralEquilibrium(testCase)
[cars,~] = carConfig();
car = cars{1,1};
P = [4,0.15,20,0.2,0.45,0,0,0.02,0.02];
[~,~,latResidual] = car.equations(P);
[~,ceq] = car.constraint9(P,P(3)*P(5),P(3)/P(5));

verifyEqual(testCase,numel(ceq),8);
verifyEqual(testCase,ceq(2),latResidual,'AbsTol',1e-10);
end

function testVectorizedTireForceMatchesScalarCorners(testCase)
car = carConfigBaseline();
car.gripScaleF = 0.91;
car.gripScaleR = 0.87;
steerLeft = 5.2;
steerRight = 5.0;
alpha = [-2.1;1.4;-3.2;2.2];
kappa = [0.01;-0.02;0.04;-0.03];
Fz = [620;580;710;680];
gamma = [-1.2;1.0;-0.8;0.7];

[Fx,Fy,Fxw] = car.tireForce(steerLeft,steerRight,alpha,kappa,Fz,gamma);
scales = [car.gripScaleF;car.gripScaleF;car.gripScaleR;car.gripScaleR];
wheelFx = zeros(4,1);
wheelFy = zeros(4,1);
for k = 1:4
    wheelFx(k) = scales(k)*car.tire.F_x(alpha(k),kappa(k),Fz(k),gamma(k));
    wheelFy(k) = scales(k)*car.tire.F_y(alpha(k),kappa(k),Fz(k),gamma(k));
end
steer = [steerLeft;steerRight;0;0];
expectedFx = wheelFx.*cosd(steer)-wheelFy.*sind(steer);
expectedFy = wheelFx.*sind(steer)+wheelFy.*cosd(steer);

verifyEqual(testCase,Fx,expectedFx,'AbsTol',1e-10);
verifyEqual(testCase,Fy,expectedFy,'AbsTol',1e-10);
verifyEqual(testCase,Fxw,wheelFx(1:2),'AbsTol',1e-10);
end

function testPitchTransferDiagnosticUsesConsistentBodyForce(testCase)
car = carConfigBaseline();
P = [0.1,0.2,10,0.1,0.3,0,0,0.01,0.01];
[~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,ssInfo] = car.equations(P);

verifyEqual(testCase,ssInfo.long_load_transfer, ...
    ssInfo.longitudinal_force_for_load_transfer*car.h_g/car.W_b, ...
    'AbsTol',1e-10);
verifyTrue(testCase,isfinite(ssInfo.longitudinal_force_for_load_transfer));
end

function [P,engineRpm] = stateBetweenConfiguredAndLegacyRedline(car)
speeds = linspace(10,40,31);
kappas = linspace(0,0.4,401);
for speed = speeds
    for kappa = kappas
        omega = (kappa+1)/car.R*speed;
        [rpm,~] = car.powertrain.engine_rpm(omega,omega,speed);
        if rpm > car.powertrain.redline + 100 && rpm < 13000
            P = [0,1,speed,0,0,repmat(kappa,1,4)];
            engineRpm = rpm;
            return
        end
    end
end
error('test_qssHardening:noRedlineProbe', ...
    'Could not construct a QSS state between configured and legacy redlines.');
end
