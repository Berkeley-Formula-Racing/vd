function [car_cell,sampleTable] = parameters_loop(cP,aP,eP,DTp,Bp,tP,samplingType,numSamples)
%PARAMETERS_LOOP Build lap and acceleration cars for grids or DOE samples.
%   FullFactorial treats every vector as explicit levels. LHS and Random
%   treat every non-scalar vector as [minimum, maximum] bounds. LHS uses a
%   base-MATLAB stratified sampler, so Statistics Toolbox is not required.

if nargin < 7 || isempty(samplingType), samplingType = "FullFactorial"; end
if nargin < 8, numSamples = []; end

names = { ...
    'mass','driver_weight','accel_driver_weight','wheelbase','weight_dist', ...
    'track_width','wheel_radius','cg_height','roll_center_height_front', ...
    'roll_center_height_rear','R_sf','I_zz','cda','cla','distribution', ...
    'cla_p_deg_p','D_p_deg_p','accel_cda','accel_cla','acc_cla_p_deg_p', ...
    'acc_D_p_deg_p','redline','shift_point','shift_time','final_drive', ...
    'drivetrain_efficiency','G_d1','G_d2_overrun','G_d2_driving', ...
    'brake_distribution','max_braking_torque','gamma_f','gamma_r','p_i', ...
    'ackermann','camber_compliance_f','camber_compliance_r','static_r_toe', ...
    'grip_scaling_front','grip_scaling_rear','I_wheel','I_driveline','Crr'};

values = { ...
    cP.mass,cP.driver_weight,cP.accel_driver_weight,cP.wheelbase,cP.weight_dist, ...
    cP.track_width,cP.wheel_radius,cP.cg_height,cP.roll_center_height_front, ...
    cP.roll_center_height_rear,cP.R_sf,cP.I_zz,aP.cda,aP.cla,aP.distribution, ...
    aP.cla_p_deg_p,aP.D_p_deg_p,aP.accel_cda,aP.accel_cla,aP.acc_cla_p_deg_p, ...
    aP.acc_D_p_deg_p,eP.redline,eP.shift_point,eP.shift_time,DTp.final_drive, ...
    DTp.drivetrain_efficiency,DTp.G_d1,DTp.G_d2_overrun,DTp.G_d2_driving, ...
    Bp.brake_distribution,Bp.max_braking_torque,tP.gamma_f,tP.gamma_r,tP.p_i, ...
    cP.ackermann,cP.camber_compliance_f,cP.camber_compliance_r,cP.static_r_toe, ...
    tP.grip_scaling_front,tP.grip_scaling_rear,cP.I_wheel,cP.I_driveline,cP.Crr};

P = sampleValues(values,samplingType,numSamples);
sampleTable = array2table(P,'VariableNames',names);
numRuns = size(P,1);
car_cell = cell(numRuns,2);

for i = 1:numRuns
    q = struct();
    for j = 1:numel(names)
        q.(names{j}) = P(i,j);
    end

    aero = Aero(q.cda,q.cla,q.distribution,q.cla_p_deg_p,q.D_p_deg_p);
    accelAero = Aero(q.accel_cda,q.accel_cla,q.distribution, ...
        q.acc_cla_p_deg_p,q.acc_D_p_deg_p);
    powertrain = Powertrain(q.redline,q.shift_point,eP.gears, ...
        eP.primary_reduction,eP.torque_fn,q.shift_time,q.final_drive, ...
        q.wheel_radius,q.drivetrain_efficiency,q.G_d1,q.G_d2_overrun, ...
        q.G_d2_driving,q.brake_distribution,q.max_braking_torque);
    tire = Tire2(q.p_i,tP.Fx_parameters,tP.Fy_parameters, ...
        tP.friction_scaling_factor);

    common = {q.wheelbase,q.weight_dist,q.track_width,q.wheel_radius, ...
        q.cg_height,q.roll_center_height_front,q.roll_center_height_rear, ...
        q.R_sf,q.I_zz,q.gamma_f,q.gamma_r,q.camber_compliance_f, ...
        q.camber_compliance_r};
    tail = {powertrain,tire,q.ackermann,q.static_r_toe,q.grip_scaling_front, ...
        q.grip_scaling_rear,q.I_wheel,q.I_driveline,q.Crr};

    car_cell{i,1} = Car(q.mass+q.driver_weight,common{:},aero,tail{:});
    car_cell{i,2} = Car(q.mass+q.accel_driver_weight,common{:},accelAero,tail{:});
end
end

function P = sampleValues(values,samplingType,numSamples)
nVar = numel(values);
switch lower(string(samplingType))
    case {"fullfactorial","grid"}
        grid = cell(1,nVar);
        [grid{:}] = ndgrid(values{:});
        P = zeros(numel(grid{1}),nVar);
        for j = 1:nVar, P(:,j) = grid{j}(:); end

    case {"lhs","random"}
        validateattributes(numSamples,{'numeric'},{'scalar','integer','positive'});
        U = rand(numSamples,nVar);
        if strcmpi(samplingType,"LHS")
            for j = 1:nVar
                U(:,j) = (randperm(numSamples)' - U(:,j))/numSamples;
            end
        end
        P = zeros(numSamples,nVar);
        for j = 1:nVar
            bounds = values{j};
            if isscalar(bounds)
                P(:,j) = bounds;
            else
                lo = min(bounds); hi = max(bounds);
                P(:,j) = lo + U(:,j)*(hi-lo);
            end
        end

    otherwise
        error('parameters_loop:badSamplingType', ...
            'samplingType must be FullFactorial, LHS, or Random.');
end
end
