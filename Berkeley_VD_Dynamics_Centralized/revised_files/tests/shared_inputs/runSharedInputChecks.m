function results = runSharedInputChecks()
%RUNSHAREDINPUTCHECKS Execute regression checks in MATLAB or GNU Octave.
% Run from any working directory after adding this file's folder to the path.
% MATLAB additionally builds both production car entrypoints. Octave tests
% the shared inputs, Tire2 force equations, Car equation samples, and a
% headless lateral-response integration. It does not run Simulink or the apps.
here = fileparts(mfilename('fullpath'));
repo = fileparts(fileparts(here));
root = fullfile(repo,'Full Car Models');
oldPath = path;
cleanup = onCleanup(@() path(oldPath)); %#ok<NASGU>
addpath(root,fullfile(root,'config'),fullfile(root,'carComponents'), ...
    fullfile(root,'utilities'),fullfile(root,'Suspension Transient Model'));
isOctave = exist('OCTAVE_VERSION','builtin') ~= 0;
if isOctave, addpath(fullfile(here,'octave_compat')); end
ref = load(fullfile(here,'upstream_defaults.mat'));
cfg = vdVehicleConfig(); tireCfg = vdTireConfig();
actual = cfg;
actual.engine.torque_fn = actual.engine.torque_fn();
actual.aero.map_path = 'aeromap_b26.csv';
for name = {'car','aero','engine','drivetrain','brakes','events'}
    assertEqual(actual.(name{1}),ref.baseline.(name{1}));
end
actualTire = tireCfg;
actualTire.model_artifact = fullfile('Magic Formula','experimental', ...
    'lc0_nd_tire','results','lc0_nd_pacejka_comparison.mat');
for name = fieldnames(ref.baseline.tire)'
    assertEqual(actualTire.(name{1}),ref.baseline.tire.(name{1}));
end
ramp = vdRampConfig(); ramp.defaultMapPath = 'aeromap_b26.csv';
assertEqual(ramp,ref.ramp);
fprintf('PASS: main and ramp default inputs match upstream fixture\n');
count = checkHistoricalDefaults(cfg,tireCfg);
fprintf('PASS: %d historical study inputs match upstream literals\n',count);
checkSharedInputEdits(root);
fprintf('PASS: temporary master mass/pressure edits propagate to ramp and presets\n');

% Exercise file validation, including a correctly shaped but invalid fit.
missing = [tempname '.mat'];
assertError(@() vdLoadTireParameters(missing,26),'vdTireConfig:missingFile');
Xbestcell = num2cell([NaN ones(1,25)]); %#ok<NASGU>
save(missing,'Xbestcell');
deleteTemp = onCleanup(@() delete(missing)); %#ok<NASGU>
assertError(@() vdLoadTireParameters(missing,26),'vdTireConfig:badParameters');
assertError(@() vdLoadTireParameters(tireCfg.Fx_file,42),'vdTireConfig:badParameters');
fprintf('PASS: missing, non-finite and wrong-length tire fits fail clearly\n');

% Use the default constructor route, so the shared camber file is exercised.
tire = Tire2(tireCfg.p_i,tireCfg.Fx_parameters,tireCfg.Fy_parameters, ...
    tireCfg.friction_scaling_factor);
p = ref.tireInputs;
fx = tire.F_x(p(:,1),p(:,2),p(:,3),p(:,4));
fy = tire.F_y(p(:,1),p(:,2),p(:,3),p(:,4));
assertEqual(fx,ref.expectedFx); assertEqual(fy,ref.expectedFy);
assert(all(fx(p(:,3)<=0)==0) && all(fy(p(:,3)<=0)==0));
% Independently check the Octave adapter against nested 1-D interpolation.
for i = 1:size(p,1)
    a = p(i,1); g = p(i,4); c = ref.camber;
    vals = [interp1(c.camber4indices,c.camber4ratio,-a,'linear','extrap'), ...
        interp1(c.camber2indices,c.camber2ratio,-a,'linear','extrap'),1, ...
        interp1(c.camber2indices,c.camber2ratio,a,'linear','extrap'), ...
        interp1(c.camber4indices,c.camber4ratio,a,'linear','extrap')];
    expected = interp1([-4 -2 0 2 4],vals,g,'linear','extrap');
    assertEqual(tire.camberInterpolator(a,g),expected);
end
fprintf('PASS: %d tire states match upstream Fx/Fy and camber interpolation\n',size(p,1));

% The acceleration car uses static aero in the production parameter loop.
% Build that same numerical core through the shared production constructor.
car = accelerationCore(cfg,tireCfg,tire);
assertEqual(car.M,cfg.car.mass+cfg.car.accel_driver_weight);
assertEqual(car.l_f+car.l_r,cfg.car.wheelbase);
for study = {'modal','decoupled','transient'}
    suspension = vdApplySuspensionConfig(car,study{1});
    assertEqual(suspension.k,cfg.suspension.(study{1}).k);
    assertEqual(suspension.Ixx,cfg.suspension.(study{1}).Ixx);
    if ~strcmp(study{1},'vibration')
        [M,C,K] = calc_MCK(suspension,zeros(1,4),zeros(1,4));
        assert(all(isfinite([M(:);C(:);K(:)])));
    end
end
% The vibration script already assigns an undeclared Car.c property upstream.
% Check its relocated inputs on a struct without changing the production class.
vibration = vdApplySuspensionConfig(struct(),'vibration');
assertEqual(vibration.c,cfg.suspension.vibration.c);
assertEqual(vibration.k,cfg.suspension.vibration.k);
fprintf('PASS: three Car suspension presets, 7-DOF matrices and vibration inputs\n');

% Fixed inputs sample the actual Car.equations path. This checks execution,
% not the engineering accuracy of the upstream equations.
samples = 0;
for speed = [5 15 25]
    for steer = [-2 0 2]*pi/180
        P = [steer .1 speed 0 steer*speed/car.W_b 0 0 .02 .02];
        [rpm,beta,ay,ax,yawacc,wheelacc] = car.equations(P);
        assert(all(isfinite([rpm beta ay ax yawacc wheelacc(:)'])));
        samples = samples+1;
    end
end
% A small imposed-speed steering response using the same Car equations.
[t,x] = ode45(@(time,state) responseRhs(time,state,car),linspace(0,2,101),[0;0]);
assert(all(isfinite(x(:))) && max(abs(x(:,1)))>0);
fprintf('PASS: %d Car equation states and %d response samples\n',samples,numel(t));

results = struct('runtime',version,'tireStates',size(p,1), ...
    'maxFxDifference',max(abs(fx-ref.expectedFx)), ...
    'maxFyDifference',max(abs(fy-ref.expectedFy)), ...
    'equationStates',samples,'responseTime',t,'responseStates',x, ...
    'MATLAB_entrypoints_executed',false);
if ~isOctave
    addpath(genpath(root));
    [cars,eventParams] = carConfig();
    [rampCar,rampConfig] = carConfigBaseline();
    assertEqual(cars{1,1}.M,cfg.car.mass+cfg.car.driver_weight);
    assertEqual(cars{1,2}.M,car.M);
    assertEqual(rampCar.M,car.M);
    assertEqual(eventParams,cfg.events);
    assertEqual(rampConfig.fixedParameters.vehicle,ramp.fixedParameters.vehicle);
    results.MATLAB_entrypoints_executed = true;
    fprintf('PASS: MATLAB production entrypoints\n');
else
    fprintf('NOT RUN: MATLAB table/aeromap entrypoints, apps and Simulink\n');
end
end

function car = accelerationCore(cfg,t,tire)
c = cfg.car; a = cfg.aero; e = cfg.engine; d = cfg.drivetrain; b = cfg.brakes;
aero = Aero(a.accel_cda,a.accel_cla,a.distribution,a.acc_cla_p_deg_p,a.acc_D_p_deg_p);
powertrain = Powertrain(e.redline,e.shift_point,e.gears,e.primary_reduction, ...
    e.torque_fn(),e.shift_time,d.final_drive,c.wheel_radius,d.drivetrain_efficiency, ...
    d.G_d1,d.G_d2_overrun,d.G_d2_driving,b.brake_distribution,b.max_braking_torque);
q = c;
for name = {'gamma_f','gamma_r','grip_scaling_front','grip_scaling_rear'}
    q.(name{1}) = t.(name{1});
end
kin = struct();
for name = {'roll_gradient_deg_per_g','rear_roll_camber_outer_deg_per_deg', ...
        'rear_roll_camber_inner_deg_per_deg','ride_camber_front_deg_per_in', ...
        'ride_camber_rear_deg_per_in'}
    kin.(name{1}) = c.(name{1});
end
kin.camberModelData = load(t.camber_model_file,'betaL','betaR');
car = buildCarFromParameterSet(q,c.accel_driver_weight,aero,powertrain,tire, ...
    [],kin,c.R_sf);
end

function dx = responseRhs(time,state,car)
speed = 10; steer = (time>=.2)*pi/180;
P = [steer 0 speed state(1) state(2) 0 0 0 0];
[~,~,ay,~,yawacc] = car.equations(P);
% Car.equations returns the body-frame lateral velocity derivative.
dx = [ay; yawacc];
end

function assertEqual(a,b)
if isstruct(b)
    assert(isequal(sort(fieldnames(a)),sort(fieldnames(b))));
    for name = fieldnames(b)', assertEqual(a.(name{1}),b.(name{1})); end
elseif isnumeric(b) || islogical(b)
    assert(isequal(size(a),size(b)));
    assert(all(abs(a(:)-b(:)) <= 1e-11*max(1,abs(b(:)))));
else
    assert(strcmp(char(a),char(b)));
end
end
function assertError(fn,id)
try
    fn();
catch err
    assert(strcmp(err.identifier,id));
    return
end
error('Expected error %s.',id);
end
