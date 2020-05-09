%% Description
% RIR_Habets.m
% An example script for the Habets RIR generator

close all;


%% Room specification
c = 340;                    % Sound velocity (m/s)
fs = 16000;                 % Sample frequency (samples/s)
r = [3.5 1.5 2];            % Receiver position [x y z] (m)
s = [2.5 3 2];              % Source position [x y z] (m)
L = [4 5 3];                % Room dimensions [x y z] (m) % one meter away from the wall
beta = 0.7;                 % Reverberation time RT60 (s)
%beta = [0.671 0.671 0.671 0.671 0.671 0.671];  % Reflection coefficients x 6
n = 4096*3;                   % Number of samples

% Generate room impulse response
% https://www.audiolabs-erlangen.de/fau/professor/habets/software/rir-generator
h = rir_generator(c, fs, r, s, L, beta, n);


%% Generate output signal

[x, Fs] = audioread('resources/IEEE_sentences/ieee01f06.wav');
y = filter(h, 1, x);        % Output with the same length as input

t_x_vals = (0:length(x)-1)/Fs;
t_h_vals = (0:n-1)/fs;
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
