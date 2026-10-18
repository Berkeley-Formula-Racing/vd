clear all;close all;clc

LF_Weight =	[48.3	57.7	68.3	83.7	90	98.3	114.7	183	196.4]-48.3;
RF_Weight =	[49.3	59.5	69	85.4	90.2	101.3	117.2	186	200.2]-49.3;
LR_Weight =	[53.3	62.9	72.5	71	76.2	83	82.1	105	135]-53.3;
RR_Weight =	[53.2	64.2	73.9	73.3	77	87.4	86.1	114	140.4]-53.2;

LFH = [0	0.002	0.004	0.006	0.008	0.01	0.017	0.031	0.034];
RFH = [0	0.002	0.004	0.007	0.008	0.01	0.011	0.029	0.033];
LRBH = [0	0.005	0.0075	0.0075	0.009	0.012	0.011	0.021	0.032];
RRBH = [0	0.004	0.007	0.007	0.0085	0.012	0.011	0.021	0.038];

x = linspace(0,0.040,1000);

figure
scatter(LFH,LF_Weight)
LF_stiffness = LFH'\LF_Weight';
hold on
plot(x,LF_stiffness*x)
xlabel('Corner Deflection (in)')
ylabel('Change in Corner Weight (lb)')
title('Front Left')
legend(['Stiffness = ' num2str(round(LF_stiffness,0)) ' (lb/in)'],'Location','Northwest')

figure
scatter(RFH,RF_Weight)
RF_stiffness = RFH'\RF_Weight';
hold on
plot(x,RF_stiffness*x)
xlabel('Corner Deflection (in)')
ylabel('Change in Corner Weight (lb)')
title('Front Right')
legend(['Stiffness = ' num2str(round(RF_stiffness,0)) ' (lb/in)'],'Location','Northwest')

figure
scatter(LRBH,LR_Weight)
LR_stiffness = LRBH'\LR_Weight';
hold on
plot(x,LR_stiffness*x)
xlabel('Corner Deflection (in)')
ylabel('Change in Corner Weight (lb)')
title('Rear Left')
legend(['Stiffness = ' num2str(round(LR_stiffness,0)) ' (lb/in)'],'Location','Northwest')

figure
scatter(RRBH,RR_Weight)
RR_stiffness = RRBH'\RR_Weight';
hold on
plot(x,RR_stiffness*x)
xlabel('Corner Deflection (in)')
ylabel('Change in Corner Weight (lb)')
title('Rear Right')
legend(['Stiffness = ' num2str(round(RR_stiffness,0)) ' (lb/in)'],'Location','Northwest')

%%


front_stiffness = round(mean([LF_stiffness RF_stiffness]),-1)
rear_stiffness = round(mean([LR_stiffness RR_stiffness]),-1)


  