%% Description
% CDR_MicSpacing.m
% Average CDR vs. Microphone Spacing (2 mics)
% Fixed config for room, source, source-receiver distance along x-axis
% Average CDR obtained across all frequency and time indices
% Two ways of changing the microphone spacing (Configuration 1 and 2)

% References:
% [1] Andreas Schwarz, Walter Kellermann, "Coherent-to-Diffuse Power Ratio Estimation for Dereverberation", 
% IEEE/ACM Trans. on Audio, Speech and Lang. Proc., 2015 (under review); preprint available: arXiv:1502.03784
% PDF: http://arxiv.org/pdf/1502.03784
% [2] Emanuël Habets, https://www.audiolabs-erlangen.de/fau/professor/habets/software/rir-generator

% Dependencies:
% 1. RIR generator
% 2. Functions from the Schwarz demo code
       
close all;

addpath(genpath('resources/Schwarz_lib'));

%% Initial parameters and configuration

L = [6.61 5.11 2.95];                                       % Room dimensions [x y z] (m)
s = [1.1 2.5 1.5];                                          % Source position [x y z] (m)
c = 342;                                                    % Speed of sound (m/s)
fs = 16000;                                                 % Sampling frequency (Hz, samples/s)
n = 12288;                                                  % Number of samples
T_60 = 1.6;                                                 % Reverberation time (s)

TxRxDistance = [0.5 2 4];                                   % Source-receiver distance (m)

% Different microphone spacings tested
spacing = [0.04 0.06 0.08 0.10 0.15 0.20 0.25 0.3 0.35 0.4 0.45 0.5 0.6 0.7 0.8 1];

% Initialize variables for results visualization
CDR_avg = zeros(length(TxRxDistance), length(spacing));
r1_all = zeros(length(TxRxDistance), length(spacing), 3);
r2_all = zeros(length(TxRxDistance), length(spacing), 3);

% Filterbank initialization
K = 512;                                                    % FFT size
N = 128;                                                    % Frame shift
Lp = 1024;                                                  % Prototype filter length
load('resources/Schwarz_lib/filterbank/prototype_K512_N128_Lp1024.mat');

% Other parameters for CDR estimation
lambda = 0.68;                                              % Forgetting factor for PSD estimation
mu = 1.3;                                                   % Noise oversubtraction factor
G_floor = 0.1;                                              % Gain floor

% Alpha and beta used for spectral subtraction
%ss_alpha = 1; ss_beta = 1;                                 % Power subtraction 
ss_alpha = 2; ss_beta = 0.5;                                % Magnitude subtraction
%ss_alpha = 2; ss_beta = 1;                                 % Wiener filter

% Chosen estimator
estimator = @estimate_cdr_nodoa;                            % DOA-independent estimator (CDRshwarz3)

% Load clean speech signal
[sig, fs_in] = audioread('resources/clean_speech/ieee01f05.wav');


%% Obtain average CDR

for i = 1:length(TxRxDistance)
    for j = 1:length(spacing)

        d_mic = spacing(j);
        
        % Configuration 1
        r1 = [s(1)+TxRxDistance(i) s(2) 1.5];               % Mic 1 position [x y z] (m)
        r2 = [s(1)+TxRxDistance(i) s(2)+d_mic 1.5];         % Mic 2 position [x y z] (m)

        % Configuration 2
        % r1 = [s(1)+TxRxDistance(i) s(2)-d_mic/2 1.5];     % Mic 1 position [x y z] (m)
        % r2 = [s(1)+TxRxDistance(i) s(2)+d_mic./2 1.5]; 	% Mic 2 position [x y z] (m)
        
        r1_all(i,j,:) = r1;
        r2_all(i,j,:) = r2;

        h1 = rir_generator(c, fs, r1, s, L, T_60, n);
        h2 = rir_generator(c, fs, r2, s, L, T_60, n);

        % Create 2-channel reverberant speech
        x1 = filter(h1, 1, sig);
        x2 = filter(h2, 1, sig);
        x_total = [x1 x2];
        x = resample(x_total, fs, fs_in);
        
        X=DFTAnaRealEntireSignal(x, K, N, p);               % Analysis filterbank

        Pxx = estimate_psd(X, lambda);
        Cxx = estimate_cpsd(X(:,:,1), X(:,:,2), lambda)./ sqrt(Pxx(:,:,1).* Pxx(:,:,2));
        frequency = linspace(0, fs/2, K/2+1)'; 
        Cnn = sinc(2 * frequency * d_mic/c);

        % Apply CDR estimator
        CDR = estimator(Cxx, Cnn);
        CDR = max(real(CDR), 0);
        CDR_avg(i,j) = mean(CDR, 'all');                    % Averaged over time and frequency indices 
        
        weights = spectral_subtraction(CDR,ss_alpha,ss_beta,mu);
        weights = max(weights, G_floor);
        weights = min(weights,1);
        Postfilter_input = sqrt(mean(abs(X).^2,3)) .* exp(1j*angle(X(:,:,1)));
        Processed = weights .* Postfilter_input;

        y = DFTSynRealEntireSignal(Processed, K, N,p);      % Synthesis filterbank
      
    end

end

%% Average CDR vs. Microphone Separation
 
figure('position',[0 0 600 450]);
scatter(spacing, 10*log10(CDR_avg(1,:)), 70, 'r', 'o', 'filled');
title('Average CDR vs. Microphone Spacing');
xlabel('Microphone Spacing/m');
ylabel('Average CDR/dB');
hold on
scatter(spacing, 10*log10(CDR_avg(2,:)), 70, 'b', '*');
scatter(spacing, 10*log10(CDR_avg(3,:)), 70, 'm', 'd', 'filled');
hold off
grid on;
legend({'Close', 'Mid-distance', 'Far'}, 'Location', 'southeast');
set(findall(gcf,'type','axes'),'fontsize',16);
set(findall(gcf,'type','text'),'fontSize',22);
% fig = gcf;
% fig.PaperPositionMode = 'auto';
% savefig(fig, 'CDR-vs-MicSpacing1.fig');
% saveas(fig, 'CDR-vs-MicSpacing1.png');  


figure('position',[0 0 600 450]);
plot3(s(:,1),s(:,2),s(:,3),'k*', ...
    r1_all(1,:,1), r1_all(1,:,2), r1_all(1,:,3), 'r+', r2_all(1,:,1), r2_all(1,:,2), r2_all(1,:,3), 'r^', ...
    r1_all(2,:,1), r1_all(2,:,2), r1_all(2,:,3), 'b+', r2_all(2,:,1), r2_all(2,:,2), r2_all(2,:,3), 'b^', ...
    r1_all(3,:,1), r1_all(3,:,2), r1_all(3,:,3), 'm+', r2_all(3,:,1), r2_all(3,:,2), r2_all(3,:,3), 'm^');
grid on;
legend({'Source', 'Mic 1 (close)', 'Mic 2 (close)', 'Mic 1 (mid)', 'Mic 2 (mid)', 'Mic 1 (far)', 'Mic 2 (far)'})
xlabel('x');
ylabel('y');
zlabel('z');
xlim([0 L(1)])
ylim([0 L(2)])
zlim([0 L(3)])
title('Positions of Source and Microphones');
set(findall(gcf,'type','axes'),'fontsize',16);
set(findall(gcf,'type','text'),'fontSize',22);
% fig = gcf;
% fig.PaperPositionMode = 'auto';
% savefig(fig, 'CDR-vs-MicSpacing-positions1.fig');
% saveas(fig, 'CDR-vs-MicSpacing-positions1.png');



