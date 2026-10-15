%% Car Parameters
clear; %close all; clc

entrypointDir = fileparts(mfilename('fullpath'));
run(fullfile(entrypointDir, 'bootstrap.m'));

% input car parameters
carCell = carConfig();
car = carCell{1,1};
car = vdApplySuspensionConfig(car,'decoupled');

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
