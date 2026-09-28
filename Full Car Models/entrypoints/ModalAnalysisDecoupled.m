%% Car Parameters
clear; %close all; clc

entrypointDir = fileparts(mfilename('fullpath'));
run(fullfile(entrypointDir, 'bootstrap.m'));

load(fullfile(fullCarModelsDir, 'Suspension Transient Model', 'dampingcurves.mat'));
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
[m,c,k] = calc_MCK_decoupled(car,tirePos,tireVel);
[V,D] = eig(inv(m)*k);

% sort eigenvalues/vectors
[d,ind] = sort(diag(D));
Ds = D(ind,ind);
Vs = V(:,ind);

% natural frequencies
omega_n = diag(sqrt(Ds))/(2*pi);

%figure
hold on
for n = 1:7
    subplot(7,1,n)
    [maximum,index] = max(abs(Vs(:,n)));
    plot(Vs(:,n)/maximum*sign(Vs(index,n)),'--o');
    title(['\omega_n =' num2str(round(omega_n(n),2)) ' \zeta =' num2str(round(zeta(n),2))])
    set(gcf,'Position',[184 66 600 714]);
    set(gca,'XTick',1:8,'XTickLabel',{'Bounce','Roll Front', 'Roll Rear','Pitch',...
        'FL Tire' 'FR Tire' 'RL Tire' 'RR Tire'})
end
