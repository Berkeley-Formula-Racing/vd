close all;clc;
setup_paths
%% Temperature and Slip Angle vs Time
%{
load('A1654raw24.mat')
figure
%ax1 = subplot(3,1,1);
plot(ET,TSTC)
hold on

load('A1965raw6.mat')
plot(ET,TSTC)

load('A1965raw15.mat')
plot(ET,TSTC)

xlim([0 500])
xlabel('Time (s)')
ylabel('Tire Center Surface Temperature (F)')
legend('18x7.5-10 R25B','16x7.5-10 R25B','16x7.5-10 LC0')

%%
ax2 = subplot(3,1,2);
plot(ET,SA)

xlim([0 90])
xlabel('Time (s)')
ylabel('Slip Angle (deg)')

ax3 = subplot(3,1,3);
load('A1654raw24.mat')
plot(ET,FY)
hold on

load('A1965raw6.mat')
plot(ET,FY)

load('A1965raw15.mat')
plot(ET,FY)
xlim([0 90])
linkaxes([ax2,ax3],'x')%linkaxes([ax1,ax2,ax3],'x'

%%
figure
load('A1654raw24.mat')
figure
plot(SA(1:10000),FY(1:10000))
hold on

load('A1965raw6.mat')
plot(SA(1:10000),FY(1:10000))

%xlim([0 300])
xlabel('Time (s)')
ylabel('Tire Center Surface Temperature (F)')
legend('18x7.5-10 R25B','16x7.5-10 R25B')

%% Good Stuff
 
figure
load('A1654raw24.mat')

plot(ET,FY)
hold on
load('A1965raw6.mat')

plot(ET,FY)
xlim([0 300])

hold on
load('A1965raw15.mat')
plot(ET,FY)

figure
plot(ET,TSTC)
xlim([0 300])

figure
plot(TSTC,FY)

%% Good Stuff 2
figure
ET_input = [0 90];
SA_input = -11;

file_name = "C:\Users\johny\Desktop\CS\vd\Magic Formula\TTC Documentation\Round 9\RawData_Cornering_Matlab_USCS_Round9_Runs1to15\A2356raw8.mat";
[SA_out, FY_out, FZ_out, P_out, TSTC_out, ET_out, index] = TireTemperatureParser2(ET_input, SA_input, file_name);
scatter(ET_out,FY_out)
hold on
%{
file_name = 'A1965raw6.mat';
[SA_out, FY_out, FZ_out, P_out, TSTC_out, ET_out, index] = TireTemperatureParser2(ET_input, SA_input, file_name);
scatter(ET_out,FY_out)
%}
file_name = 'A1965raw15.mat';
[SA_out, FY_out, FZ_out, P_out, TSTC_out, ET_out, index] = TireTemperatureParser2(ET_input, SA_input, file_name);
scatter(ET_out,FY_out)

xlabel('Time (s)')
ylabel('Lateral Force (lb)')
legend('16x7.5-10 R20','16x7.5-10 LC0')
%}
%% Good Stuff 3
ET_input = [0 200]; % Time range
SA_inputs = 4:1:12; % Slip angles to iterate over

% File paths for datasets
file_name_1 = "C:\Users\johny\Desktop\CS\vd\Magic Formula\TTC Documentation\Round 9\RawData_Cornering_Matlab_USCS_Round9_Runs1to15\A2356raw8.mat";
file_name_2 = 'A1965raw15.mat';

for i = 1:length(SA_inputs)
    SA_input = SA_inputs(i);

    % Parse data for first dataset
    [SA_out, FY_out_1, FZ_out_1, P_out_1, TSTC_out_1, ET_out_1, IA_out_1, NFY_out_1, index_1] = TireTemperatureParser2(ET_input, SA_input, file_name_1);

    % Parse data for second dataset
    [SA_out, FY_out_2, FZ_out_2, P_out_2, TSTC_out_2, ET_out_2, NFY_out_2, IA_out_2, index_2] = TireTemperatureParser2(ET_input, SA_input, file_name_2);

    % Create a new figure for each slip angle
    figure;
    hold on;

    % Plot first dataset (R20) in red
    scatter(TSTC_out_1, smooth(-FY_out_1, 1), 'r', 'DisplayName', '16x7.5-10 R20');

    % Plot second dataset (LC0) in blue
    scatter(TSTC_out_2, smooth(-FY_out_2, 1), 'b', 'DisplayName', '16x7.5-10 LC0');

    % Customize plot
    xlabel('Temperature (F)');
    ylabel('Lateral Force (lbs)');
    title(sprintf('Lateral Force vs Temperature (SA = %d)', SA_input));
    legend('show');
    grid on;
    hold off;
    filename = sprintf('R20-LC0-Temp-SA%d.png', SA_input);
    saveas(gcf, filename);
end


%{
%% Good Stuff 4

figure
ET_input = [77 80];
SA_input = 0;

file_name = 'A1654raw24.mat';
[SA_out, FY_out, FZ_out, P_out, TSTC_out, ET_out, IA_out, index] = TireTemperatureParser2(ET_input, SA_input, file_name);
scatter(SA_out,FY_out)
hold on

ET_input = [40 43];
SA_input = 0;
file_name = 'A1965raw6.mat';
[SA_out, FY_out, FZ_out, P_out, TSTC_out, ET_out, IA_out, index] = TireTemperatureParser2(ET_input, SA_input, file_name);
scatter(SA_out,FY_out)

ET_input = [30 33];
SA_input = 0;

file_name = 'A1965raw15.mat';
[SA_out, FY_out, FZ_out, P_out, TSTC_out, ET_out,  IA_out,index] = TireTemperatureParser2(ET_input, SA_input, file_name);
scatter(SA_out,FY_out)
xlabel('Slip Angle (deg)')
ylabel('Lateral Force (Lb)')
legend('18x7.5-10 R25B','16x7.5-10 R25B','16x7.5-10 LC0')

figure
ET_input = [74 77];
SA_input = 0;

file_name = 'A1654raw24.mat';
[SA_out, FY_out, FZ_out, P_out, TSTC_out, ET_out,  IA_out,index] = TireTemperatureParser2(ET_input, SA_input, file_name);
plot(TSTC_out)
hold on
save('TemperatureCompensated_18R25B.mat','SA_out', 'FY_out', 'FZ_out', 'P_out', 'TSTC_out', 'ET_out',  'IA_out','index')

ET_input = [40 43];
SA_input = 0;
file_name = 'A1965raw6.mat';
[SA_out, FY_out, FZ_out, P_out, TSTC_out, ET_out,  IA_out,index] = TireTemperatureParser2(ET_input, SA_input, file_name);
plot(TSTC_out)
save('TemperatureCompensated_16R25B.mat','SA_out', 'FY_out', 'FZ_out', 'P_out', 'TSTC_out', 'ET_out',  'IA_out','index')

ET_input = [30 33];
SA_input = 0;

file_name = 'A1965raw15.mat';
[SA_out, FY_out, FZ_out, P_out, TSTC_out, ET_out, IA_out, index] = TireTemperatureParser2(ET_input, SA_input, file_name);
plot(TSTC_out)
xlabel('Index')
ylabel('Temperature (F)')
legend('18x7.5-10 R25B','16x7.5-10 R25B','16x7.5-10 LC0')
save('TemperatureCompensated_16LC0.mat','SA_out', 'FY_out', 'FZ_out', 'P_out', 'TSTC_out', 'ET_out', 'IA_out', 'index')


%%
hold on
load('A1965raw6.mat')

plot(ET,FY)
xlim([0 90])

hold on
load('A1965raw15.mat')
plot(ET,FY)

figure
plot(ET,TSTC)
xlim([0 90])

figure
plot(TSTC,FY)

%%
load('A1654run24.mat')
figure
subplot(2,1,1)
plot(ET,TSTC)
hold on

load('A1965run6.mat')
plot(ET,TSTC)

%xlim([0 300])
xlabel('Time (s)')
ylabel('Tire Center Surface Temperature (F)')
legend('18x7.5-10 R25B','16x7.5-10 R25B')

subplot(2,1,2)
%plot(ET,SA)

xlim([0 300])
xlabel('Time (s)')
ylabel('Slip Angle (deg)')

%% Good Stuff 4
P_input = 12;
IA_input = 0;
FZ_input = [250];
SA = 10;
file_name = 'A1654run24.mat';

[SA_out, FY_out, FZ_out, P_out, TSTC_out, ET_out, index] = TireTemperatureParser(P_input, IA_input, FZ_input, SA, file_name);

plot(SA_out,TSTC_out)
hold on
file_name = 'A1965run15.mat';

[SA_out, FY_out, FZ_out, P_out, TSTC_out, ET_out, index] = TireTemperatureParser(P_input, IA_input, FZ_input, SA, file_name);
plot(SA_out,TSTC_out)
xlabel('Slip Angle (deg)')
ylabel('Tire Center Surface Temp (F)')
legend('18x7.5-10 R25B','16x7.5-10 LC0')
title('FZ = 250, P = 12',...
    'FontSize',15);



%%
file_name = '"C:\Users\johny\Desktop\CS\vd\Magic Formula\TTC Documentation\Round 9\RawData_Cornering_Matlab_USCS_Round9_Runs1to15 (2)\A2356raw8.mat"';
[SA_out, FY_out, FZ_out, P_out, TSTC_out, ET_out, index] = TireTemperatureParser(P_input, IA_input, FZ_input, SA, file_name);

plot(SA_out,FY_out)
hold on
file_name = 'A1965run6.mat';

[SA_out, FY_out, FZ_out, P_out, TSTC_out, ET_out, index] = TireTemperatureParser(P_input, IA_input, FZ_input, SA, file_name);
plot(SA_out,FY_out)

%%
figure
subplot(2,1,1)
plot(ET_out,TSTC_out)
subplot(2,1,2)
plot(ET_out,SA_out)

%%
P_input = 10;
IA_input = 0;
FZ_input = 150;
SA = 12;

[SA_out, FY_out, FZ_out, P_out, TSTC_out, ET_out, index] = TireTemperatureParser(P_input, IA_input, FZ_input, SA, file_name);
plot(SA_out,TSTC_out)

figure
plot(TSTC_out,-FY_out)

%%
figure
ET_input = [0 10000];
SA_input = [-11];
P_input = 12;
IA_input = 0;
FZ_input = [50 150 250];

file_name = "C:\Users\johny\Desktop\CS\vd\Magic Formula\TTC Documentation\Round 9\RawData_Cornering_Matlab_USCS_Round9_Runs1to15\A2356raw8.mat";
[SA_out, FY_out, FZ_out, P_out, TSTC_out, ET_out, index] = TireTemperatureParser3(ET_input,P_input, IA_input, FZ_input, SA_input, file_name);

scatter(TSTC_out,FY_out);
hold on


file_name = 'A1965raw15.mat';
[SA_out, FY_out, FZ_out, P_out, TSTC_out, ET_out, index] = TireTemperatureParser3(ET_input,P_input, IA_input, FZ_input, SA_input, file_name);

scatter(TSTC_out,FY_out);
xlabel('Temperature (F)')
ylabel('Friction Coefficient')
legend('16x7.5-10 R20','16x7.5-10 LC0')

%}
