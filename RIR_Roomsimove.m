%% Description
% RIR_Roomsimove.m
% An example script for the Roomsimove RIR generator
% All configuration details (e.g. sourcer, sensor, sampling frequency etc.)
% need to be specified in the 'mic_config.txt' file

close all;
addpath(genpath('resources/helper_scripts'));

%% Generate outout signal
% Parameter specification is stored in a text file
% http://homepages.loria.fr/evincent/software/Roomsimove.zip

[x, Fs] = audioread('resources/clean_speech/ieee01f05.wav');
h = roomsimove_single('mic_config.txt', [2.5; 3; 2]);   % Sensor config, Source config
y = fftfilt(h, x);
h_fs = 16000;                                           % From the config text file
t_x_vals = (0:length(x)-1)/Fs;
t_h_vals = (0:n-1)/fs;
[h_edc] = edc(h);

%% Visualisation

figure('position',[0 0 600 450]);
plot(t_x_vals, x);
xlabel('Time/s');
ylabel('Amplitude');
grid on;
title('Clean Speech Signal');
set(findall(gcf,'type','axes'),'fontsize',16);
set(findall(gcf,'type','text'),'fontSize',22);
% fig = gcf;
% fig.PaperPositionMode = 'auto';
% savefig(fig, 'Figures/RIR/Habets-clean.fig');
% saveas(fig, 'Figures/RIR/Habets-clean.png');

figure('position',[0 0 600 450]);
plot(t_h_vals, h);
xlabel('Time/s');
ylabel('Amplitude');
grid on;
title('Room Impulse Response');
set(findall(gcf,'type','axes'),'fontsize',16);
set(findall(gcf,'type','text'),'fontSize',22);
% fig = gcf;
% fig.PaperPositionMode = 'auto';
% savefig(fig, 'Figures/RIR/Habets-h.fig');
% saveas(fig, 'Figures/RIR/Habets-h.png');

figure('position',[0 0 600 450]);
plot(t_h_vals, h_edc);
set(gca,'FontSize', 14);
title("Energy Decay Curve of the Room Impulse Response", 'Fontsize', 18);
xlabel('Time/s', 'Fontsize', 20);
ylabel('Power (dB)', 'Fontsize', 20);
grid on;

figure('position',[0 0 600 450]);
plot(t_x_vals, y);
xlabel('Time/s');
ylabel('Amplitude');
grid on;
title('Reverberant Signal');
set(findall(gcf,'type','axes'),'fontsize',16);
set(findall(gcf,'type','text'),'fontSize',22);
% fig = gcf;
% fig.PaperPositionMode = 'auto';
% savefig(fig, 'Figures/RIR/Habets-reverberant.fig');
% saveas(fig, 'Figures/RIR/Habets-reverberant.png');

%soundsc(y, Fs);
%audiowrite('out.wav',y,Fs);