%% Description
% trueCDR.m
% Calculate the ground truth of CDR

% CDR depends on the complex spatial coherence
% Complex spatial coherence uses auto-power and cross-power spectra

close all;

%% Generate signals received at the microphones

% Room specification
c = 342;            % Speed of sound (m/s)
fs = 16000;         % Sampling frequency (samples/s)
L = [5 5 3];        % Room dimensions [x y z] (m)
s = [2.5 3 1.5];    % Source position [x y z] (m)
% 2 mics
d_mic = 0.06;       % Mic spacing
r = [3 2 1.5; 3 2+d_mic 1.5];   % Receiver position [x y z] (m)
T60 = 1.2;          % Reverberation time (s)
n = 12288;          % Number of samples

% RIR simulations using the Habets RIR generator
% https://www.audiolabs-erlangen.de/fau/professor/habets/software/rir-generator
h = rir_generator(c, fs, r, s, L, T60, n);

% Load clean speech as input signal
[s_in, fs_in] = audioread('resources/speech/ieee01f05.wav');

% Convolve with RIR to obtain the reverberant signal (mixed field)
%x = filter(h, 1, s_in);   % only one channel at a time
x = fftfilt(h', s_in);
x = resample(x, fs, fs_in);

load('lib/filterbank/prototype_K512_N128_Lp1024.mat');

%% Complex spatial coherence

% Mixed-signal (i.e. signal received at the mics)
X = DFTAnaRealEntireSignal(x, 512, 128, p);
Cxx = cpsd(X(:,:,1), X(:,:,2)) ./ cpsd(X(:,:,1), X(:,:,1));

% Direct-path component
% Keep all other config the same except now T60 = 0
h_d = rir_generator(c, fs, r, s, L, 0, n);
x_d = fftfilt(h_d', s_in);
x_d = resample(x_d, fs, fs_in);
X_d = DFTAnaRealEntireSignal(x_d, 512, 128, p);
Css = cpsd(X_d(:,:,1), X_d(:,:,2)) ./ cpsd(X_d(:,:,1), X_d(:,:,1));

% Diffuse noise component (reverberation)
% method 1: assuming signal and noise mutually orthogonal 
% auto- and cross-psd
% apsd_n = cpsd(X(:,:,1), X(:,:,1)) - cpsd(X_d(:,:,1), X_d(:,:,1));
% cpsd_n = cpsd(X(:,:,1), X(:,:,2)) - cpsd(X_d(:,:,1), X_d(:,:,2));
% Cnn = cpsd_n ./ apsd_n;
% method 2: use h
h_n = h - h_d;
x_n = fftfilt(h_n', s_in);
x_n = resample(x_n, fs, fs_in);
X_n = DFTAnaRealEntireSignal(x_n, 512, 128, p);
Cnn = cpsd(X_n(:,:,1), X_n(:,:,2)) ./ cpsd(X_n(:,:,1), X_n(:,:,1));

%% CDR calculation
CDR = (Cnn - Cxx) ./ (Cxx - Css);
CDR_avg = mean(CDR, 'all');
CDR_avg_db = 10*log10(CDR_avg);
