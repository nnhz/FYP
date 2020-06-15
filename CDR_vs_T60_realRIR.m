%% Description
% CDR_vs_Rooms.m
% To investigate the relationship between CDR and T60 using real RIRs taken
% from the ACE Challenge Corpus
% CDR estimation code taken from [1]

% Reference:
% [1] Andreas Schwarz, Walter Kellermann, "Coherent-to-Diffuse Power Ratio
% Estimation for Dereverberation", IEEE/ACM Trans. on Audio, Speech and
% Lang. Proc., 2015 (under review); preprint available: arXiv:1502.03784
% PDF: http://arxiv.org/pdf/1502.03784

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

CDR_avg = zeros(1, n_rooms);
c = 342;                    
fs = 16000;                          
n = 8192; 

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
[x, fs_in] = audioread('resources/clearn_speech/ieee01f05.wav');

%% Obtain average CDR for the same speech signal in different rooms
% Using demo_cdr_dereverb.m by Andreas Schwarz (schwarz@lnt.de)
% Adapted for literature review experiment
% Estimator prop 3 (since it is DOA-independent) 

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
    
    % analysis filterbank
    X=DFTAnaRealEntireSignal(y, K, N, p);

    % estimate PSD and coherence
    Pxx = estimate_psd(X, lambda);
    Cxx = estimate_cpsd(X(:,:,1), X(:,:,2), lambda)./ sqrt(Pxx(:,:,1).* Pxx(:,:,2));

    frequency = linspace(0, fs/2, K/2+1)'; % frequency axis

    % define coherence models
    Cnn = sinc(2 * frequency * d_mic/c); % diffuse noise coherence; not required for estimate_cdr_nodiffuse

    % apply CDR estimator (=SNR)
    CDR = estimator(Cxx, Cnn);
    CDR = max(real(CDR), 0);
    CDR_avg(i) = mean(CDR, 'all');  % averaged over time and frequency indices 

    weights = spectral_subtraction(CDR,ss_alpha,ss_beta,mu);
    weights = max(weights, G_floor);
    weights = min(weights,1);

    % postfilter input is computed from averaged PSDs of both microphones
    Postfilter_input = sqrt(mean(abs(X).^2,3)) .* exp(1j*angle(X(:,:,1)));

    % apply postfilter
    Processed = weights .* Postfilter_input;

    % synthesis filterbank
    z_cdr = DFTSynRealEntireSignal(Processed, K, N,p);
    
end

%% Average CDR vs. T60 (with real RIRs)

figure('position',[0 0 600 450]);
scatter(T_60, 10*log10(CDR_avg), 80, 'filled');
title('Average CDR vs. T60');
xlabel('T60/s');
ylabel('Average CDR/dB');
grid on;
set(findall(gcf,'type','axes'),'fontsize',16);
set(findall(gcf,'type','text'),'fontSize',22);
fig = gcf;
fig.PaperPositionMode = 'auto';
% savefig(fig, 'CDR-vs-T60.fig');
% saveas(fig, 'CDR-vs-T60.png');



