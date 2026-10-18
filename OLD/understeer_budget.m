% K % C coefficient conventions are (+) for all understeer effects and (-) for oversteer effects regardless of % front or rear designation.

% Force & Moment based K & C coefficients are 'per 1000 N' or 'per 100 Nm' to avoid the ridiculously small
% values likely to be obfuscated by those having exceptional brilliance.

% Unit conventions are N, Nm, deg, kg.
%----------------------------------------------
% Cornering compliance budget after simulation has completed.
clear

CAF = 136*4.448 ; CAR = 136*4.448% Tire cornering stiffnesses per tire (N/deg)
CCF = 0.313 ; CCR = 0.323 % Tire camber stiffnesses per tire (N/deg)
NAF = 50.082 ; NAR = 42.897 % Tire aligning moment stiffnesses per tire (Nm/deg)
CGF = 47.657 ; CGR = 43.392 % Tire camber stiffnesses per tire (N/deg)
NGF = 4.766 ; NGR = 4.339 % Tire alinging moment camber stiffnesses per tire (Nm/deg)

m = 168;
g = 9.8;


L = 62 * 25.4 % WHEELBASE (mm)
weight_dist = 0.51; % rearwards

a = L*weight_dist; % front axle to cg, m
b = L*(1-weight_dist); % rear axle to cg, m

FYFperg = m * g * b / L% front axle sideforce per g from simulation (N/g)
FYRperg = 8430.9 % Rear
NFperg = -403.16 % Front axle aligning moment per g from simulation (Nm/g)
NRperg = -317.93 % Rear

WF = m*(1-weight_dist); WR = 885.000 % AXLE MASS (kg)
WUF = 17.07; WUR = 99.000 %UNSPRUNG MASS (kg)

EF = 0.001; ER = 0.055 %ROLL STEER (deg steer/deg roll)
GF = 0.789; GR = -0.720 %ROLL CAMBER (deg camber/deg roll)
KF = 2288; KR = 367; %ROLL STIFFNESS (front roll stiffness (Nm/deg)

EYF = 0.191; EYR = 0.045 %L.F. STEER (deg/1000 N)
GYF = 0.291; GYR = -0.328 %L.F. CAMBER (deg/1000 N)
ENF = 0.774; ENR = -0.113 %A.T. STEER (deg/100 Nm)
GNF = 0.022; GNR = 0.023 %A.T. CAMBER (deg/100 Nm)

cc_WGT_f = 9.806*WF/2/CAF; % Cornering Compliances from weight contributions
cc_WGT_r = -9.806*WR/2/CAR;

cc_RBAT_f = -(NFperg + NRperg)/2/L/CAF; % Cornering Compliances from rigid body Mz (F)
cc_RBAT_r = -(NFperg + NRperg)/2/L/CAR; % Rear

ROLLpg = 4.47; % Roll gradient (deg/g)

cc_RS_f = ROLLpg * EF ; % Front roll steer contribution (deg/g)
cc_RS_r = ROLLpg * ER ; % Rear

cc_RC_f = ROLLpg * GF * CGF/CAF; % Roll Camber Front (deg/g)
cc_RC_r = ROLLpg * GR * CGR/CAR; % Rear

cc_LFDS_f = EYF *(FYFperg-9.806* WUF )/2/1000; % From Front L.F. deflection steer
cc_LFDS_r = EYR *(FYRperg-9.806* WUR )/2/1000; % Rear

cc_LFDC_f = GYF *CGF*(FYFperg-9.806* WUF )/2/1000/CAF ; % Front L.F. Deflection Camber
cc_LFDC_r = GYR *CGR*(FYRperg-9.806* WUR )/2/1000/CAR ; % Rear

cc_ATDS_f = - ENF *NFperg/2/100; % Front Aligning Moment Deflection Steer
cc_ATDS_r = - ENR *NRperg/2/100; % Rear

cc_ATDC_f = - GNF *CGF*NFperg/2/CAF/100; % Front Aligning Moment Deflection Camber
cc_ATDC_r = - GNR *CGR*NRperg/2/CAR/100; % Rear


DF = cc_WGT_f + cc_RBAT_f + cc_RS_f + cc_RC_f + cc_LFDS_f + cc_LFDC_f + cc_ATDS_f + cc_ATDC_f; % Front Cornering Compliance Sum
DR = cc_WGT_r + cc_RBAT_r + cc_RS_r + cc_RC_r + cc_LFDS_r + cc_LFDC_r + cc_ATDS_r + cc_ATDC_r; % Rear Cornering Compliance Sum

K = DF + DR; % DR is a negative number (always oversteering). Total Vehicle Understeer (deg/g)

cc(1,:)= [cc_WGT_f cc_RBAT_f cc_RS_f cc_RC_f cc_LFDS_f cc_LFDC_f cc_ATDS_f cc_ATDC_f];
cc(2,:)= [cc_WGT_r cc_RBAT_r cc_RS_r cc_RC_r cc_LFDS_r cc_LFDC_r cc_ATDS_r cc_ATDC_r];
% Bar Plot of Cornering Compliance Breakdown
figure('numberTitle','Off','MenuBar','None','Name','Cornering Compliance Element Summary');

% Bar chart
h = bar(cc, 'stacked');
colormap('hsv')

% Legend and labels
L = legend(h, {'Wgt/Tires','Rigid Body Mz','Roll Steer','Roll Camber', ...
    'Fy Steer','Fy Camber','Mz Steer','Mz Camber'});
set(L,'FontSize',8,'Location','North');
legend('boxoff')

% Axis setup
ylabel('Cornering Compliance Values (deg/g)')
set(gca, 'XTickLabel', {'Front', 'Rear'})
ylim([-1.5 1.5])  % Adjust based on actual values
grid on

% Reference lines
yline(1.5, '--k', 'Max Ref');
yline(0, '-k');

% Bottom text
text(1, -1.35, 'Front', 'HorizontalAlignment', 'center')
text(2, -1.35, 'Rear', 'HorizontalAlignment', 'center')
text(1.5, -1.45, 'Computed From 1^o Linear Tire Stiffnesses', 'HorizontalAlignment', 'center', 'FontSize', 8)

% Title with summary
title(['DF = ' num2str(round(DF,3)) '  |  DR = ' ...
       num2str(round(DR,3)) '  |  K = ' ...
       num2str(round(K,3)) ' deg/g'])
