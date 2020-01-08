%% Room specification
c = 340;                    % Sound velocity (m/s)
fs = 48000;                 % Sampling frequency (samples/s)
L = [3.22 5.1 2.94];        % Room dimensions [x y z] (m)
s = [1.41 1.73 1.19];       % Source position [x y z] (m)
beta = 0.4;                 % Reverberation time (s)
n = 4096;                   % Number of samples

%% Generate RIR for an ULA (with different receiver position r
% https://www.audiolabs-erlangen.de/fau/professor/habets/software/rir-generator

n_microphones = 2;                  % Number of microphones in the ULA
spacing = 0.02;                     % Space between two adjcent microphones
h = zeros(n_microphones, n);        % RIRs
r = zeros(n_microphones);

for i = 0:n_microphones-1
    r_i = [1.25 2.81+i*spacing 1.19];   % Receiver position [x y z] (m)
    h(i+1,:) = rir_generator(c, fs, r_i, s, L, beta, n);
end

figure;
plot(h(1,:));
title("Microphone 1 RIR");

figure;
plot(h(2,:));
title("Microphone 2 RIR");

%% Import the actual recorded RIR (ACE_Corpus)
[recorded, recorded_fs] = audioread('resources/ACE_Corpus/Single_803_1_RIR.wav');

figure;
plot(recorded);
title("Actual RIR");

%% Convolution to get output signals
[x, Fs] = audioread('resources/IEEE_sentences/ieee01f06.wav');

y1 = filter(h(1,:)', 1, x);     % Output same length as input
y2 = filter(h(2,:)', 1, x);
y_recorded = filter(recorded, 1, x);

%% Results

% soundsc(y1, Fs);
% soundsc(y2, Fs);
% soundsc(y_recorded, Fs);
