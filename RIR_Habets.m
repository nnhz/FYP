%% Description
% RIR_Habets.m
% An example script for the Habets RIR generator

close all;
addpath(genpath('resources/helper_scripts'));

%% Room specification
c = 340;                    % Sound velocity (m/s)
fs = 16000;                 % Sample frequency (samples/s)
r = [3.5 1.5 2];            % Receiver position [x y z] (m)
s = [2.5 3 2];              % Source position [x y z] (m)
L = [4 5 3];                % Room dimensions [x y z] (m) % one meter away from the wall
T_60 = 0.4;                 % Reverberation time (s)
% beta = [0.671 0.671 0.671 0.671 0.671 0.671];  % Reflectivity
n = 12288;                   % Number of samples

% Generate room impulse response
% https://www.audiolabs-erlangen.de/fau/professor/habets/software/rir-generator
h = rir_generator(c, fs, r, s, L, T_60, n);

%% Generate output signal

[x, Fs] = audioread('resources/clean_speech/ieee01f05.wav');
y = filter(h, 1, x);        
t_x_vals = (0:length(x)-1)/Fs;
t_h_vals = (0:n-1)/fs;
[h_edc] = edc(h);

%% Visualization of Results

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
