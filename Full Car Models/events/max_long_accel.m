function [x_accel,long_accel,long_accel_guess,diagnostics] = ...
    max_long_accel(long_vel_guess,car,x0,solverOptions)
% uses fmincon to minimize the objective function subject to constraints
% optimizes longitudinal acceleration with given lateral acceleration constraint
% disp('max long accel');
if nargin < 3 || isempty(x0) % no initial guess supplied
    %initial guesses
    steer_angle_guess = 0;
    throttle_guess = 1;
    lat_vel_guess = 0;
    yaw_rate_guess = 0;

    kappa_1_guess = 0;
    kappa_2_guess = 0;
    kappa_3_guess = 0.01;
    kappa_4_guess = 0.01;
    
    x0 = [steer_angle_guess,throttle_guess,long_vel_guess,lat_vel_guess,yaw_rate_guess,kappa_1_guess,...
        kappa_2_guess,kappa_3_guess,kappa_4_guess];
end
if nargin < 4 || isempty(solverOptions)
    solverOptions = struct();
end
if ~isstruct(solverOptions) || ~isscalar(solverOptions)
    error('max_long_accel:invalidSolverOptions', ...
        'solverOptions must be a scalar struct.');
end

x0(3) = long_vel_guess;

% bounds
steer_angle_bounds = [0,0];
throttle_bounds = [0,1]; 
long_vel_bounds = [long_vel_guess,long_vel_guess];
lat_vel_bounds = [0,0];
yaw_rate_bounds = [0,0];
kappa_1_bounds = [0,0];
kappa_2_bounds = [0,0];
kappa_3_bounds = [0,0.2];
kappa_4_bounds = [0,0.2];

A = [];
b = [];
Aeq = [1 0 0 0 0 0 0 0 0
       0 0 1 0 0 0 0 0 0
       0 0 0 1 0 0 0 0 0
       0 0 0 0 1 0 0 0 0
       0 0 0 0 0 1 0 0 0
       0 0 0 0 0 0 1 0 0
       0 0 0 0 0 0 0 1 -1]; % symmetric rear slip: kappa_4 = kappa_3
beq = [0 long_vel_guess 0 0 0 0 0];
lb = [steer_angle_bounds(1),throttle_bounds(1),long_vel_bounds(1),lat_vel_bounds(1),...
    yaw_rate_bounds(1),kappa_1_bounds(1),kappa_2_bounds(1),kappa_3_bounds(1),kappa_4_bounds(1)];
ub = [steer_angle_bounds(2),throttle_bounds(2),long_vel_bounds(2),lat_vel_bounds(2),...
    yaw_rate_bounds(2),kappa_1_bounds(2),kappa_2_bounds(2),kappa_3_bounds(2),kappa_4_bounds(2)];

% objective function: longitudinal acceleration (forwards)
f = @(P) -car.long_accel(P);

% no lateral acceleration constraint
constraint = @(P) car.constraint1(P);

% default algorithm is interior-point
maxFunctionEvaluations = getSolverOption(solverOptions, ...
    'maxFunctionEvaluations',2000);
constraintTolerance = getSolverOption(solverOptions, ...
    'constraintTolerance',1e-2);
stepTolerance = getSolverOption(solverOptions,'stepTolerance',1e-10);
display = getSolverOption(solverOptions,'display','off');
if isstring(display)
    display = char(display);
end
options = optimoptions('fmincon', ...
    'MaxFunctionEvaluations',maxFunctionEvaluations, ...
    'ConstraintTolerance',constraintTolerance, ...
    'StepTolerance',stepTolerance,'Display',display);

% fval: objective function value (v^2/r)
% exitflag meaning: 1 = converged, 2 = change in x less than step tolerance
%   (optimality condition not fulfilled, but solution still found
%   0 = function evaluations exceeded (not converging)
%   -2 = no feasible point found 
[x,fval,exitflag,~,~,~,~] = fmincon(f,x0,A,b,Aeq,beq,lb,ub,constraint,options);
x(9) = x(8);

long_accel_guess = x;

[engine_rpm,beta,~,long_accel,~,~,omega,current_gear,...
~,Fz,alpha,T] = car.equations(x);

% generate table of control variable values
x_accel = [exitflag long_accel x(3)*x(5) x omega(1:4) engine_rpm current_gear beta...
    Fz(1:4) alpha(1:4) T(1:4)];
[~] = generate_table(x_accel);

long_accel = -fval;

if nargout >= 4
    diagnostics = struct();
    diagnostics.state = x;
    diagnostics.exitflag = exitflag;
    [diagnostics.c,diagnostics.ceq] = car.constraint1(x);
    diagnostics.max_inequality_violation = max([diagnostics.c(:);0]);
    diagnostics.max_equality_residual = max(abs(diagnostics.ceq(:)));
    diagnostics.metrics = car.metrics(x);
    diagnostics.pure_ay0 = abs(diagnostics.metrics.gLat) <= 1e-12;
end
end

function value = getSolverOption(solverOptions,name,default)
value = default;
if isfield(solverOptions,name) && ~isempty(solverOptions.(name))
    value = solverOptions.(name);
end
end
