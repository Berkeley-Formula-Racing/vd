% K % C coefficient conventions are (+) for all understeer effects and (-) for oversteer effects regardless of % front or rear designation.

% Force & Moment based K & C coefficients are 'per 1000 N' or 'per 100 Nm' to avoid the ridiculously small
% values likely to be obfuscated by those having exceptional brilliance.

% Unit conventions are N, Nm, deg, kg.
%----------------------------------------------
% Cornering compliance budget after simulation has completed.
clear

addpath(fullfile(fileparts(mfilename('fullpath')),'Full Car Models','config'));
sharedStudyVehicle = vdVehicleConfig();
sharedStudyTire = vdTireConfig([],[],false);
CAF = sharedStudyTire.studies.understeer_budget.CAF; CAR = sharedStudyTire.studies.understeer_budget.CAR% Tire cornering stiffnesses per tire (N/deg)
CCF = sharedStudyTire.studies.understeer_budget.CCF; CCR = sharedStudyTire.studies.understeer_budget.CCR % Tire camber stiffnesses per tire (N/deg)
NAF = sharedStudyTire.studies.understeer_budget.NAF; NAR = sharedStudyTire.studies.understeer_budget.NAR % Tire aligning moment stiffnesses per tire (Nm/deg)
CGF = sharedStudyTire.studies.understeer_budget.CGF; CGR = sharedStudyTire.studies.understeer_budget.CGR % Tire camber stiffnesses per tire (N/deg)
NGF = sharedStudyTire.studies.understeer_budget.NGF; NGR = sharedStudyTire.studies.understeer_budget.NGR % Tire alinging moment camber stiffnesses per tire (Nm/deg)

m = sharedStudyVehicle.studies.understeer_budget.m;
g = 9.8;


L = sharedStudyVehicle.studies.understeer_budget.L% WHEELBASE (mm)
weight_dist = sharedStudyVehicle.studies.understeer_budget.weight_dist; % rearwards

a = L*weight_dist; % front axle to cg, m
b = L*(1-weight_dist); % rear axle to cg, m

FYFperg = m * g * b / L% front axle sideforce per g from simulation (N/g)
FYRperg = sharedStudyVehicle.studies.understeer_budget.FYRperg % Rear
NFperg = sharedStudyVehicle.studies.understeer_budget.NFperg % Front axle aligning moment per g from simulation (Nm/g)
NRperg = sharedStudyVehicle.studies.understeer_budget.NRperg % Rear

WF = m*(1-weight_dist); WR = sharedStudyVehicle.studies.understeer_budget.WR % AXLE MASS (kg)
WUF = sharedStudyVehicle.studies.understeer_budget.WUF; WUR = sharedStudyVehicle.studies.understeer_budget.WUR %UNSPRUNG MASS (kg)

EF = sharedStudyVehicle.studies.understeer_budget.EF; ER = sharedStudyVehicle.studies.understeer_budget.ER %ROLL STEER (deg steer/deg roll)
GF = sharedStudyVehicle.studies.understeer_budget.GF; GR = sharedStudyVehicle.studies.understeer_budget.GR %ROLL CAMBER (deg camber/deg roll)
KF = sharedStudyVehicle.studies.understeer_budget.KF; KR = sharedStudyVehicle.studies.understeer_budget.KR; %ROLL STIFFNESS (front roll stiffness (Nm/deg)

EYF = sharedStudyVehicle.studies.understeer_budget.EYF; EYR = sharedStudyVehicle.studies.understeer_budget.EYR %L.F. STEER (deg/1000 N)
GYF = sharedStudyVehicle.studies.understeer_budget.GYF; GYR = sharedStudyVehicle.studies.understeer_budget.GYR %L.F. CAMBER (deg/1000 N)
ENF = sharedStudyVehicle.studies.understeer_budget.ENF; ENR = sharedStudyVehicle.studies.understeer_budget.ENR %A.T. STEER (deg/100 Nm)
GNF = sharedStudyVehicle.studies.understeer_budget.GNF; GNR = sharedStudyVehicle.studies.understeer_budget.GNR %A.T. CAMBER (deg/100 Nm)

cc_WGT_f = 9.806*WF/2/CAF; % Cornering Compliances from weight contributions
cc_WGT_r = -9.806*WR/2/CAR;

cc_RBAT_f = -(NFperg + NRperg)/2/L/CAF; % Cornering Compliances from rigid body Mz (F)
cc_RBAT_r = -(NFperg + NRperg)/2/L/CAR; % Rear

ROLLpg = sharedStudyVehicle.studies.understeer_budget.ROLLpg; % Roll gradient (deg/g)

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
