%% Description
% RIR_Roomsimove.m
% An example script for the Roomsimove RIR generator
% All configuration details (e.g. sourcer, sensor, sampling frequency etc.)
% need to be specified in the 'mic_config.txt' file

close all;


%% Generate outout signal
% Parameter specification is stored in a text file
% http://homepages.loria.fr/evincent/software/Roomsimove.zip

[x, Fs] = audioread('resources/IEEE_sentences/ieee01f06.wav');
h = roomsimove_single('mic_config.txt', [2.5; 3; 2]);  % Sensor config, Source config
y = fftfilt(h, x);

% TODO: rewrite to not have hardcoded value
h_fs = 16000;  % From the config text file

t_x_vals = (0:length(x)-1)/Fs;
t_h_vals = (0:length(h)-1)/h_fs;
h_db = mag2db(h);
h_60db = max(h_db) - 60;


%% Visualisation

figure;
plot(t_x_vals, x);
set(gca,'FontSize', 12);
title("Input Signal: Clean Speech", 'Fontsize', 20);
xlabel('Time/s', 'Fontsize', 20);
ylabel('Amplitude', 'Fontsize', 20);

figure;
plot(t_h_vals, h);
set(gca,'FontSize', 12);
title("Room Impulse Response", 'Fontsize', 20);
xlabel('Time/s', 'Fontsize', 20);
ylabel('Amplitude', 'Fontsize', 20);

figure;
semilogy(t_h_vals, h_db);
yline(h_60db, '--', '-60dB', 'Color', 'r');
set(gca,'FontSize', 12);
title("Room Impulse Response", 'Fontsize', 20);
xlabel('Time/s', 'Fontsize', 20);
ylabel('Magnitude (dB)', 'Fontsize', 20);

figure;
plot(t_x_vals, y);
set(gca,'FontSize', 12);
title("Output Signal: Noisy Speech", 'Fontsize', 20);
xlabel('Time/s', 'Fontsize', 20);
ylabel('Amplitude', 'Fontsize', 20);

%soundsc(y, Fs);
%audiowrite('out.wav',y,Fs);
