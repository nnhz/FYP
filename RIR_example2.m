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
t_h_vals = (0:length(H)-1)/h_fs;

figure;
plot(t_x_vals, x);
title("Input Signal: Clean Speech");
xlabel('Time/s');
ylabel('Amplitude');

figure;
plot(t_h_vals, H);
title("Room Impulse Response");
xlabel('Time/s');
ylabel('Amplitude');

figure;
plot(t_x_vals, y);
title("Output Signal: Noisy Speech");
xlabel('Time/s');
ylabel('Amplitude');