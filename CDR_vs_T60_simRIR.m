%% Description
% CDR_vs_T60_simRIR.m
% To investigate the relationship between CDR and T60 using generated RIRs
% Positions of source and microphones, room size is fixed
% Only T60 is changed

close all;

addpath(genpath('resources/Schwarz_lib'));

%% Initial parameters and configuration
d_mic = 0.06;
L = [6.61 5.11 2.95];       % Room dimensions [x y z] (m)
% L = [10 7 3];       

s = [3 2.5 1.5];          % Source position [x y z] (m) T_60 = 0.332;
r1 = [7 3.48 1.5];
r2 = [7 3.48+d_mic 1.5];

% s = [5 3.5 1.5];          
% r1 = [1.1 2.5 1.5];
% r2 = [1.1 2.5+d_mic 1.5];

vol = prod(L);

n_rooms = 19;
T_60 = linspace(0.2, 2, n_rooms);
% n_rooms = 39;
% T_60 = linspace(0.2, 4, n_rooms);

Dc = zeros(1, n_rooms);

% Obtain critical distance Dc
for i=1:n_rooms
    Dc(i) = 0.1 * sqrt( (1*vol) / (pi * T_60(i)));
end

Dc_bound = max(Dc);

c = 342;                    % Speed of sound (m/s)
fs = 16000;                 % Sampling frequency (Hz, samples/s)
n = 4096*3;                 % Number of samples

% Initialize variables for visualisation of results
DRR_all = zeros(n_rooms, 2);
CDR_avg = zeros(1, n_rooms);


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
[x, fs_in] = audioread('resources/clean_speech/ieee01f05.wav');

%% Obtain average DRR
% Using estDRRFromRIR from the ACE Challenge Corpus
% DRR can be estimated using 1 mic, but for comparison with CDR, 2 mics are used (and take the average)

for i = 1:n_rooms

    % Generate RIRs, assuming 2 mirophones
    % https://www.audiolabs-erlangen.de/fau/professor/habets/software/rir-generator

    h1 = rir_generator(c, fs, r1, s, L, T_60(i), n);
    h2 = rir_generator(c, fs, r2, s, L, T_60(i), n);

    DRR_all(i, :) = estDRRFromRIR([h1;h2]', fs);

    % Create 2-channel noisy speech
    y1 = filter(h1, 1, x);
    y2 = filter(h2, 1, x);
    y_total = [y1 y2];

    y = resample(y_total, fs, fs_in);

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
    CDR_avg(i) = mean(CDR, 'all'); % average over time and frequency indices 
    
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



%% Results

DRR_avg = mean(DRR_all, 2);

figure('position',[0 0 600 450]);
scatter(T_60, DRR_avg, 80, 'filled');
title('Average DRR vs. T60');
xlabel('T60/s');
ylabel('Average DRR/dB');
grid on;
set(findall(gcf,'type','axes'),'fontsize',16);
set(findall(gcf,'type','text'),'fontSize',22);
fig = gcf;
fig.PaperPositionMode = 'auto';
% savefig(fig, 'DRR-vs-T60.fig');
% saveas(fig, 'DRR-vs-T60.png');

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


figure('position',[0 0 600 450]);
plot3(s(:,1), s(:,2), s(:,3), 'k*', r1(:,1), r1(:,2), r1(:,3), 'b^', r2(:,1), r2(:,2), r2(:,3), 'r>');
grid on;
legend({'Source', 'Mic 1', 'Mic 2'})
xlim([0 L(1)])
ylim([0 L(2)])
zlim([0 L(3)])
title('Positions of Source and Microphones');
set(findall(gcf,'type','axes'),'fontsize',16);
set(findall(gcf,'type','text'),'fontSize',22);
fig = gcf;
fig.PaperPositionMode = 'auto';
% savefig(fig, 'CDRDRR-vs-T60-positions.fig');
% saveas(fig, 'CDRDRR-vs-T60-positions.png');

% audiowrite('y_example.wav', y, fs);
% audiowrite('z_example.wav', z, fs);

