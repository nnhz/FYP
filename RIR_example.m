close all;

%% Room specification
c = 340;                    % Sound velocity (m/s)
fs = 16000;                 % Sample frequency (samples/s)
r = [3.5 1.5 2];              % Receiver position [x y z] (m)
s = [2.5 3 2];              % Source position [x y z] (m)
L = [4 5 3];              % Room dimensions [x y z] (m) % one meter away from the wall
beta = 0.4;                 % Reverberation time (s)
n = 4096;                   % Number of samples

% Generate room impulse response
% https://www.audiolabs-erlangen.de/fau/professor/habets/software/rir-generator
h = rir_generator(c, fs, r, s, L, beta, n);

%% Generate output signal
[x, Fs] = audioread('resources/IEEE_sentences/ieee01f06.wav');
y = filter(h, 1, x);        % Output with the same length as input

t_x_vals = (0:length(x)-1)/Fs;
t_h_vals = (0:n-1)/fs;

figure;
plot(t_x_vals, x);
title("Input Signal: Clean Speech");
xlabel('Time/s');
ylabel('Amplitude');

figure;
plot(t_h_vals, h);
title("Room Impulse Response");
xlabel('Time/s');
ylabel('Amplitude');

figure;
plot(t_x_vals, y);
title("Output Signal: Noisy Speech");
xlabel('Time/s');
ylabel('Amplitude');

soundsc(y, Fs);
