close all;

addpath(genpath('resources/Schwarz_lib'));


%% Initial parameters and configuration

% All dimensions taken drom ACE Challenge Corpus
n_rooms = 7;

dimensions = zeros(n_rooms, 3);
dimensions(1,:) = [3.32 4.83 2.95];
dimensions(2,:) = [3.22 5.1 2.94];
dimensions(3,:) = [6.61 5.11 2.95];
dimensions(4,:) = [10.3 9.07 2.63];
dimensions(5,:) = [6.93 9.73 3];
dimensions(6,:) = [13.6 9.29 2.94];
dimensions(7,:) = [4.47 5.13 3.18];

T_60 = [0.332; 0.39; 0.437; 0.371; 0.638; 1.22; 0.646];
vol = [47.3 48.3 99.6 246 202 370 72.9];
Dc = zeros(1, n_rooms);

s_positions = zeros(n_rooms, 3);
s_positions(1,:) = [1.1 2.2 1.5];
s_positions(2,:) = [1.1 2.2 1.5];
s_positions(3,:) = [1.1 2.5 1.5];
s_positions(4,:) = [3.0 4.5 1.5];
s_positions(5,:) = [2.0 4.5 1.5];
s_positions(6,:) = [2.0 4.5 1.5];
s_positions(7,:) = [1.1 2.2 1.5];

r1_positions = zeros(n_rooms, 3);
r1_positions(1,:) = [2.9 2.5 1.5];
r1_positions(2,:) = [2.9 3.9 1.5];
r1_positions(3,:) = [5.0 4.0 1.5];
r1_positions(4,:) = [8.0 6.2 1.5];
r1_positions(5,:) = [5.0 5.0 1.5];
r1_positions(6,:) = [11.0 5.0 1.5];
r1_positions(7,:) = [3.8 2.5 1.5];

r2_positions = [r1_positions(:,1) r1_positions(:,2)+0.02 r1_positions(:,3)];

distances = zeros(1, n_rooms);

% Obtain critical distance Dc
for i=1:n_rooms
    Dc(i) = 0.1 * sqrt( (1*vol(i)) / (pi * T_60(i)));
    distances(i) = norm( s_positions(i,:) - (r1_positions(i,:) + r2_positions(i,:))/2 );
end

SNR_avg = zeros(1, n_rooms);
c = 340;                    
fs = 16000;                          
n = 8192; 


%% Load clean speech signal
[x, fs_in] = audioread('resources/IEEE_sentences/ieee01f06.wav');


%% Obtain average CDR for the same speech signal in different rooms
% Using demo_cdr_dereverb.m by Andreas Schwarz (schwarz@lnt.de)
% Adapted for literature review experiment
% Estimator prop 3 (since it is DOA-independent) 



% filterbank initialization
cfg.K = 512; % FFT size
cfg.N = 128; % frame shift
cfg.Lp = 1024; % prototype filter length
load('resources/Schwarz_lib/filterbank/prototype_K512_N128_Lp1024.mat');

% algorithm and scenario configuration
cfg.fs = 16000;      % sampling rate [Hz]
cfg.c = 342;         % speed of sound [m/s]
cfg.d_mic = 0.08;   % mic spacing [m]

cfg.nr.lambda = 0.68; % smoothing factor for PSD estimation
cfg.nr.mu = 1.3;     % noise overestimation factor
cfg.nr.floor = 0.1;  % minimum gain
%cfg.nr.alpha = 1; cfg.nr.beta = 1; % power subtraction
cfg.nr.alpha = 2; cfg.nr.beta = 0.5; % magnitude subtraction
%cfg.nr.alpha = 2; cfg.nr.beta = 1; % Wiener filter

cfg.estimator = @estimate_cdr_nodoa;              % DOA-independent estimator (CDRprop3)


for i = 1:n_rooms
    
    % Generate RIRs, assuming 2 mirophones
    % https://www.audiolabs-erlangen.de/fau/professor/habets/software/rir-generator

    h1 = rir_generator(c, fs, r1_positions(i,:), s_positions(i,:), dimensions(i,:), T_60(i), n);
    h2 = rir_generator(c, fs, r2_positions(i,:), s_positions(i,:), dimensions(i,:), T_60(i), n);

    % Create 2-channel noisy speech
    y1 = filter(h1, 1, x);
    y2 = filter(h2, 1, x);
    y_total = [y1 y2];
    
    y = resample(y_total, cfg.fs, fs_in);
    
    % signal processing
    fprintf('Performing signal enhancement... ');tic;

    % analysis filterbank
    X=DFTAnaRealEntireSignal(y,cfg.K,cfg.N,p);

    % estimate PSD and coherence
    Pxx = estimate_psd(X,cfg.nr.lambda);
    Cxx = estimate_cpsd(X(:,:,1),X(:,:,2),cfg.nr.lambda)./sqrt(Pxx(:,:,1).*Pxx(:,:,2));

    frequency = linspace(0,cfg.fs/2,cfg.K/2+1)'; % frequency axis

    % define coherence models
    Cnn = sinc(2 * frequency * cfg.d_mic/cfg.c); % diffuse noise coherence; not required for estimate_cdr_nodiffuse

    % apply CDR estimator (=SNR)
    SNR = cfg.estimator(Cxx, Cnn);
    SNR = max(real(SNR),0);
    SNR_avg(i) = mean(SNR, 'all'); % average over time and frequency indices 
    
    weights = spectral_subtraction(SNR,cfg.nr.alpha,cfg.nr.beta,cfg.nr.mu);
    weights = max(weights,cfg.nr.floor);
    weights = min(weights,1);

    % postfilter input is computed from averaged PSDs of both microphones
    Postfilter_input = sqrt(mean(abs(X).^2,3)) .* exp(1j*angle(X(:,:,1)));

    % apply postfilter
    Processed = weights .* Postfilter_input;

    % synthesis filterbank
    z = DFTSynRealEntireSignal(Processed,cfg.K,cfg.N,p);
    fprintf('done (%.2fs).\n', toc);
    
    %audiowrite('resources/out.wav',z,cfg.fs);

    % Plot results
%     figure;
%     imagesc(10*log10(SNR))
%     set(gca,'YDir','normal')
%     caxis([-15 15])
%     colorbar
%     title(join(['Estimated CDR (=SNR) [dB] Room ', sprintf('%02d', i)]))
%     xlabel('frame index')
%     ylabel('subband index')
end

%% Average CDR vs. Room dimensions

figure;
scatter(T_60, 10*log10(SNR_avg));
title('Average CDR vs. T60');
xlabel('T60');
ylabel('Average CDR');
   
