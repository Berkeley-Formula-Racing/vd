clear; clc
% run this from the same directory as setup_paths
try 
    setup_paths
catch
    error('run this from same directory as setup_paths. Should be one directory up');
end

% input car parameters
carCell = carConfig();
car = carCell{1,1};
car = vdApplySuspensionConfig(car,'vibration');
car.TSmpc = .003; %has to be multiple of TSdyn
car.TSdyn = .0005;
car.Jm = 0; car.Jw = 1;
n = 8000; % number of timesteps


k_f = car.k;
k_r = car.k;
tireConfig = vdTireConfig([],[],false);
k_tf = car.k*tireConfig.vertical_stiffness.vibration_spring_multiplier;
k_tr = k_tf;
c_f = car.c;
c_r = car.c;
m = car.M;
m_f = m/10;
m_r = m_f;
k_rf = car.k_rf;
k_rr = car.k_rr;
Ix = car.Ixx;
Iy = car.Iyy;
a_1 = car.l_f;
a_2 = car.l_r;
b_1 = car.t_f/2;
b_2 = car.t_f/2;
w = car.t_f;

% state vector: [x phi theta x1 x2 x3 x4]'

m = diag([m Ix Iy m_f m_f m_r m_r]);

c11 = 2*c_f+2*c_r;
c21 = b_1*c_f-b_2*c_f-b_1*c_r+b_2*c_r;
c12 = c21;
c31 = 2*a_2*c_r-2*a_1*c_f;
c13 = c31;
c22 = b_1^2*c_f+b_2^2*c_f+b_1^2*c_r+b_2^2*c_r;
c32 = a_1*b_2*c_f-a_1*b_1*c_f-a_2*b_1*c_r+a_2*b_2*c_r;
c23 = c32;
c33 = 2*c_f*a_1^2+2*c_r*a_2^2;

c = [c11 c12 c13 -c_f -c_f -c_r -c_r;
    c21 c22 c23 -b_1*c_f b_2*c_f b_1*c_r -b_2*c_r;
    c31 c32 c33 a_1*c_f a_1*c_f -a_2*c_r -a_2*c_r;
    -c_f -b_1*c_f a_1*c_f c_f 0 0 0;
    -c_f b_2*c_f a_1*c_f 0 c_f 0 0;
    -c_r b_1*c_r -a_2*c_r 0 0 c_r 0;
    -c_r -b_2*c_r -a_2*c_r 0 0 0 c_r];

k11 = 2*k_f+2*k_r;
k21 = b_1*k_f-b_2*k_f-b_1*k_r+b_2*k_r;
k12 = k21;
k31 = 2*a_2*k_r - 2*a_1*k_f;
k13 = k31;
k22 = k_rf+b_1^2*k_f+b_2^2*k_f+b_1^2*k_r+b_2^2*k_r;
k32 = a_1*b_2*k_f-a_1*b_2*k_f-a_2*b_1*k_r+a_2*b_2*k_r;
k23 = k32;
k42 = -b_1*k_f-k_rf/w;
k24 = k42;
k52 = b_2*k_f+k_rf/w;
k25 = k52;
k33 = 2*k_f*a_1^2+2*k_r+a_2^2;
k44 = k_f+k_tf+k_rf/w^2;
k55 = k44;

k = [k11 k12 k13 -k_f -k_f -k_r -k_r;
    k21 k22 k23 k24 k25 b_1*k_r -b_2*k_r;
    k31 k32 k33 a_1*k_f a_1*k_f -a_2*k_r -a_2*k_r;
    -k_f k42 a_1*k_f k44 -k_rf/w^2 0 0;
    -k_f k52 a_1*k_f -k_rf/w^2 k55 0 0;
    -k_r b_1*k_r -a_2*k_r 0 0 k_r+k_tr 0;
    -k_r -b_2*k_r -a_2*k_r 0 0 0 k_r+k_tr];

F = [-car.M*9.81 500 0 0 0 0 0]';
    
x0 =[0 0 0 0 0 0 0]';
xd0 = [0 0 0 0 0 0 0 ]';

%xdd = inv(m)*(F-c*xdot-k*x);

% y = [x;xd]
% yd = [xd xdd]

[t,y] = ode45(@(t,y) dydt(t,y,m,F,c,k),[0 10],[x0; xd0]);

x = y(:,1:7);
xd = y(:,8:end);
% state vector: [x phi theta x1 x2 x3 x4]'
body = x(:,1);
phi = x(:,2);
theta = x(:,3);

function [yd, blah] = dydt(t,y,m,F,c,k) 
    x = y(1:7,1);
    xd = y(8:end,1);
    xdd = inv(m)*(F-c*xd-k*x);
    yd = [xd; xdd]; 
end



