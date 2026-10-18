%% Load Data

clear all;clc

% import mat file from motec
load('Sine1Period.mat');

% input variable names of channels 
data = {Engine_RPM,G_Force_Lat,G_Force_Long,GP_Volts_5,GP_Volts_6,...
    GP_Volts_7,GP_Volts_10,Ground_Speed_Left,Ground_Speed_Right...
    Gyro_Yaw_Velocity,Steering_Angle};
%% Inputs

% select data channel
x = data{10};

% plotting
figure1 = true; % comparison of filtered and unfiltered data
figure2 = true; % filter magnitude response
figure3 = true; % delay
figure4 = true; % magnitude response of filtered and unfiltered data

Fp = 2; % frequency at start of pass band
Fst = 3; % frequency at end of stop band
Ap = 0.001; % amount of ripple allowed in pass band (dB)
Ast = 50; % attenuation of stop band (dB)

%% Calculations

% sampling frequency
Fs = 1/(x.Time(2)-x.Time(1));
L = numel(x.Value);
f = Fs*(1:(L/2))/L;
fourier_transform = fft(x.Value);

% filter data
d = fdesign.lowpass('Fp,Fst,Ap,Ast',Fp,Fst,Ap,Ast,Fs);
Hd = design(d,'equiripple'); 

data_filtered2 = movmean(x.Value,50);
data_filtered = filter(Hd,x.Value);

% shift data by delay caused by filtering
grp = grpdelay(Hd);
delay = ceil(mean(grp));
L2 = numel(data_filtered);
fourier_transform2 = fft(data_filtered);
data_filtered(1:delay) = [];

%% Plotting

if figure1 == true
    figure
    plot(x.Time(1:end-delay),data_filtered,'LineWidth',2)
    hold on
    plot(x.Time,data_filtered2)
end

if figure2 == true
    fvtool(Hd);
end

if figure3 == true
    figure
    plot(grp)
end

if figure4 == true
    figure
    plot(f,abs(fourier_transform(1:L/2)));
    hold on
    plot(f,abs(fourier_transform2(1:L2/2)));
end
