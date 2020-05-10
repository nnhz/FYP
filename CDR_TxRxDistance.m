%% Description
% CDR_TxRxDistance.m
% Average CDR vs. Source-Receiver Distance
% Fixed config for room, source
% Vary the receiver positions (using 2 mics)
% Average CDR obtained across all frequency and time indices

close all;

addpath(genpath('resources/Schwarz_lib'));

%% Initial parameters and configuration

L = [6.61 5.11 2.95];       % Room dimensions [x y z] (m)
s = [1.1 2.5 1.5];          % Source position [x y z] (m) T_60 = 0.332;

c = 342;                    % Speed of sound (m/s)
fs = 16000;                 % Sampling frequency (Hz, samples/s)
n = 4096*3;                 % Number of samples
T_60 = 0.45;                % Reverberation time (s)
d_mic = 0.06;               % Mic spacing (m)

% Number of iterations in the x- and y-direction for 2 mic positions
x_itr = 6;
y_itr = 6;

% Initialize variables for visualisation of results
SNR_avg = zeros(1, x_itr * y_itr);
%SNR_max = zeros(1, x_itr * y_itr);
distances = zeros(1, x_itr * y_itr);
r1_all = zeros(x_itr * y_itr, 3);
r2_all = zeros(x_itr * y_itr, 3);

% Filterbank initialization
K = 512;                    % FFT size
N = 128;                    % Frame shift
%Lp = 1024;                  % Prototype filter length
load('resources/Schwarz_lib/filterbank/prototype_K512_N128_Lp1024.mat');

% Other parameters for CDR estimation
lambda = 0.68;              % Forgetting/smoothing factor for PSD estimation
mu = 1.3;                   % Noise oversubtraction/overestimation factor
G_floor = 0.1;              % Minimum gain

% Alpha and beta used for spectral subtraction
%ss_alpha = 1; ss_beta = 1;     % Power subtraction 
ss_alpha = 2; ss_beta = 0.5;    % Magnitude subtraction
%ss_alpha = 2; ss_beta = 1;     % Wiener filter

% Chosen estimator
estimator = @estimate_cdr_nodoa;   % DOA-independent estimator (CDRprop3)

%% Load clean speech signal
[x, fs_in] = audioread('resources/IEEE_sentences/ieee01f06.wav');


%% Obtain average CDR
% Using demo_cdr_dereverb.m by Andreas Schwarz (schwarz@lnt.de)
% Adapted for literature review experiment
% Estimator prop 3 (since it is DOA-independent)


for i = 1:y_itr
    for j = 1:x_itr

        % Generate RIRs, assuming 2 mirophones
        % https://www.audiolabs-erlangen.de/fau/professor/habets/software/rir-generator

        r1 = [1.5+j*0.5 1+i*0.5 1.5];          % Receiver position [x y z] (m)
        r2 = [1.5+j*0.5 1+d_mic+i*0.5 1.5];    % Receiver position [x y z] (m)

        r1_all(j+6*(i-1),:) = r1;
        r2_all(j+6*(i-1),:) = r2;

        distances(j+6*(i-1)) = norm(s - (r1 + r2)/2);

        h1 = rir_generator(c, fs, r1, s, L, T_60, n);
        h2 = rir_generator(c, fs, r2, s, L, T_60, n);

        % Create 2-channel noisy speech
        y1 = filter(h1, 1, x);
        y2 = filter(h2, 1, x);
        y_total = [y1 y2];

        y = resample(y_total, fs, fs_in);

        % signal processing
        fprintf('Performing signal enhancement... '); tic;

        % analysis filterbank
        X=DFTAnaRealEntireSignal(y, K, N, p);

        % estimate PSD and coherence
        Pxx = estimate_psd(X, lambda);
        Cxx = estimate_cpsd(X(:,:,1), X(:,:,2), lambda)./ sqrt(Pxx(:,:,1).* Pxx(:,:,2));

        frequency = linspace(0, fs/2, K/2+1)'; % frequency axis

        % define coherence models
        Cnn = sinc(2 * frequency * d_mic/c); % diffuse noise coherence; not required for estimate_cdr_nodiffuse

        % apply CDR estimator (=SNR)
        SNR = estimator(Cxx, Cnn);
        SNR = max(real(SNR), 0);
        SNR_avg(j+6*(i-1)) = mean(SNR, 'all'); % average over time and frequency indices 
        %SNR_max(j+6*(i-1)) = max(SNR, [], 'all');
        
        weights = spectral_subtraction(SNR,ss_alpha,ss_beta,mu);
        weights = max(weights, G_floor);
        weights = min(weights,1);

        % postfilter input is computed from averaged PSDs of both microphones
        Postfilter_input = sqrt(mean(abs(X).^2,3)) .* exp(1j*angle(X(:,:,1)));

        % apply postfilter
        Processed = weights .* Postfilter_input;

        % synthesis filterbank
        z = DFTSynRealEntireSignal(Processed, K, N,p);
        fprintf('done (%.2fs).\n', toc);

        %audiowrite('resources/out.wav',z,fs);

        % Plot results
    %     figure;
    %     imagesc(10*log10(SNR))
    %     set(gca,'YDir','normal')
    %     caxis([-15 15])
    %     colorbar
    %     title(join(['Estimated CDR (=SNR) [dB] with distance of ', sprintf('%02d', i*0.5)]))
    %     xlabel('frame index')
    %     ylabel('subband index')


    end

end

%% Average CDR vs. Source-Receiver Distance

figure;
scatter(distances, 10*log10(SNR_avg));
title('Average CDR vs. Source-Receiver Distance');
xlabel('Distance');
ylabel('Average CDR');
   
figure;
plot3(s(:,1), s(:,2), s(:,3), 'k*', r1_all(:,1), r1_all(:,2), r1_all(:,3), 'b^', r2_all(:,1), r2_all(:,2), r2_all(:,3), 'r>');
grid on;
legend({'s', 'r1', 'r2'})
xlim([0 L(1)])
ylim([0 L(2)])
zlim([0 L(3)])
title('Positions of Source and Microphones');

