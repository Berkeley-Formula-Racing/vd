%% Car Parameters
clear; close all; clc

setup_paths

load('dampingcurves.mat');
% input car parameters
carCell = carConfig();
car = carCell{1,1};
car.k = 250*4.45*39.37; % N/m
car.k_tf = car.k*2; % tire stiffness (N/m)
car.k_tr = car.k_tf;

% LSC-HSC, LSR_HSR
% LS ticks from fully clockwise
% HS turns from fully counter clockwise
% C12/R12 valving
% 1:5 HS sweep (0-x 0-x), x = 4.3 3 2 1 0
% 6:11 LS sweep (x-4.3 x-4.3), x = 2 4 6 10 15 25
car.c_compression = cdamp{6}([1 2 10:10:end],:); % damping curves ([in/s lbf])
car.c_rebound = rdamp{6}([1 2 10:10:end],:);

car.Ixx = 60;
car.Iyy = 82;
car.k_rf = 3000*6.474; % Nm/rad
car.k_rr = 7000*6.474; % Nm/rad
car.MR_F = MR18Front;
car.MR_R = MR18Front;

%% Modal Analysis (Undamped)
tirePos = zeros(1,4);
tireVel = zeros(1,4);
[m,c,k] = calc_MCK(car,tirePos,tireVel);
[V,D] = eig(inv(m)*k);

% sort eigenvalues/vectors
[d,ind] = sort(diag(D));
Ds = D(ind,ind);
Vs = V(:,ind);

% natural frequencies
omega_n = diag(sqrt(Ds))/(2*pi);

%% Damping Analysis
% simplified 3DOF model used (bounce,roll,pitch) for approximating damping
% ratio
% (7 DOF model has non-proportional damping so it can't be decoupled)
tirePos = zeros(1,4);
tireVel = zeros(1,4);
[m,c,k] = calc_MCK(car,tirePos,tireVel);
m = m(1:3,1:3);
c = c(1:3,1:3);
k = k(1:3,1:3);
[V,D] = eig(inv(m)*k);

% One technique for calculating damping ratios
% "A SIMPLIFIED APPROACH TO CALCULATE THE
% DAMPING RATIO MATRIX FOR MULTIPLE DEGREE OF
% FREEDOM VIBRATION SYSTEMS"
c_cr = 2*m^(0.5)*(m^(-0.5)*k*m^(-0.5))^(0.5)*m^(0.5);
eig(c-c_cr);
zeta2 = c_cr^(-0.5)*c*c_cr^(-0.5);
zeta2 = eig(zeta2);

% create modal damping matrices (decoupled)
% c_bar is diagonal, indicating proportional damping
V_norm = V*diag(sqrt(diag(V'*(m*V)).^(-1)));
c_bar = V_norm'*c*V_norm;
k_bar = V_norm'*k*V_norm;
zeta = 0.5*diag(c_bar)./diag(sqrt(k_bar));

zeta = [zeta; zeros(4,1)];
omega_d = omega_n.*(sqrt(abs(ones(size(zeta))-zeta.^2)));

%figure
hold on
for n = 1:7
    subplot(7,1,n)
    [maximum,index] = max(abs(Vs(:,n)));
    plot(Vs(:,n)/maximum*sign(Vs(index,n)),'--o');
    title(['\omega_n =' num2str(round(omega_n(n),2)) ' \omega_d ='...
        num2str(round(omega_d(n),2)) ' \zeta =' num2str(round(zeta(n),2))], 'Interpreter','tex')
    set(gcf,'Position',[184 66 600 714]);
    set(gca,'XTick',1:7,'XTickLabel',{'Bounce','Roll','Pitch',...
        'FL Tire' 'FR Tire' 'RL Tire' 'RR Tire'})
end

%% Quarter Car Model (No Tire Stiffness)
MR = interp1(car.MR_F(:,1),car.MR_F(:,2),0);
wheel_rate = car.k*MR^2;
omega_n_QC = sqrt(wheel_rate/(car.M/4))/(2*pi); %Hz
c_cr_QC = 2*sqrt(wheel_rate*(car.M/4)); %Ns/m
c_QC = (car.c_compression(2,2)-car.c_compression(1,2))/(car.c_compression(2,1)-car.c_compression(1,1))*175.126835*MR^2;
zeta_QC = c_QC/c_cr_QC;

%% Damping Ratio vs Damper Settings
for i = 1:6
    
    % LSC-HSC, LSR_HSR
    % LS ticks from fully clockwise
    % HS turns from fully counter clockwise
    % C12/R12 valving
    % 1:5 HS sweep (0-x 0-x), x = 4.3 3 2 1 0
    % 6:11 LS sweep (x-4.3 x-4.3), x = 2 4 6 10 15 25
    car.c_compression = cdamp{5+i}([1 2 10:10:end],:); % damping curves ([in/s lbf])
    car.c_rebound = rdamp{5+i}([1 2 10:10:end],:);
    
    tirePos = zeros(1,4);
    tireVel = zeros(1,4);
    [m,c,k] = calc_MCK(car,tirePos,tireVel);
    m = m(1:3,1:3);
    c = c(1:3,1:3);
    k = k(1:3,1:3);
    [V,D] = eig(inv(m)*k);
    
    % One technique for calculating damping ratios
    % "A SIMPLIFIED APPROACH TO CALCULATE THE
    % DAMPING RATIO MATRIX FOR MULTIPLE DEGREE OF
    % FREEDOM VIBRATION SYSTEMS"
    c_cr = 2*m^(0.5)*(m^(-0.5)*k*m^(-0.5))^(0.5)*m^(0.5);
    eig(c-c_cr);
    zeta2 = c_cr^(-0.5)*c*c_cr^(-0.5);
    zeta2 = eig(zeta2);
    
    % create modal damping matrices (decoupled)
    % c_bar is diagonal, indicating proportional damping
    V_norm = V*diag(sqrt(diag(V'*(m*V)).^(-1)));
    c_bar = V_norm'*c*V_norm;
    k_bar = V_norm'*k*V_norm;
    zeta = 0.5*diag(c_bar)./diag(sqrt(k_bar));
    
    zeta_vec{i} = zeta;
end

for j = 1:6
    bounce(j) = zeta_vec{j}(1);
    roll(j) = zeta_vec{j}(2);
    pitch(j) = zeta_vec{j}(3);
end

figure
plot(bounce,'-o')
hold on
plot(roll,'-o')
hold on
plot(pitch,'-o')
ylabel('Damping Ratio')
xlabel('Low Speed Damper Settings')
legend('Bounce','Pitch','Roll')

