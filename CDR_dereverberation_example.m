%% Description
% An example script that integrate the RIR generation and the
% dereverberation system using CDR 

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
addpath(genpath('resources/helper_scripts'));

%% Configuration 
d_mic = 0.08;                                                   % Mic spacing (m)
L = [6.61 5.11 2.95];                                           % Room dimensions [x y z] (m)
c = 342;                                                        % Speed of sound (m/s)
fs = 16000;                                                     % Sampling frequency (Hz, samples/s)
n = 4096*3;                                                     % Number of samples
T_60 = 0.7;                                                     % Reverberation time (s)

% s = [1.1 2.5 1.5];                                            % Source position [x y z] (m)
% r1 = [3 2.4 1.5];                                             % Receiver position [x y z] (m)
% r2 = [3 2.4+d_mic 1.5];                                       % Receiver position [x y z] (m)
s = [2.5 2.5 1.5];                                              % Source position [x y z] (m)
r1 = [4 3 1.5];                                                 % Receiver position [x y z] (m)
r2 = [4 3+d_mic 1.5];                                           % Receiver position [x y z] (m)
      
% Filterbank initialization
K = 512;                                                        % FFT size
N = 128;                                                        % Frame shift
Lp = 1024;                                                      % Prototype filter length
load('resources/Schwarz_lib/filterbank/prototype_K512_N128_Lp1024.mat');

% Other parameters for CDR estimation
lambda = 0.68;                                                  % Forgetting factor for PSD estimation
mu = 1.3;                                                       % Noise oversubtraction factor
G_floor = 0.1;                                                  % Gain floor

% Alpha and beta used for spectral subtraction
%ss_alpha = 1; ss_beta = 1;                                     % Power subtraction 
ss_alpha = 2; ss_beta = 0.5;                                    % Magnitude subtraction
%ss_alpha = 2; ss_beta = 1;                                     % Wiener filter

% Chosen estimator
estimator = @estimate_cdr_nodoa;   % DOA-independent estimator (CDRprop3)

% Load clean speech signal
[sig, fs_in] = audioread('resources/clean_speech/ieee01f05.wav');
% [sig, fs_in] = audioread('resources/roomC-2m-75deg.wav');     % Demo reverberant input

%% Obtain average CDR

h1 = rir_generator(c, fs, r1, s, L, T_60, n);
h2 = rir_generator(c, fs, r2, s, L, T_60, n);

% Create 2-channel noisy speech
x1 = filter(h1, 1, sig);
x2 = filter(h2, 1, sig);
x_total = [x1 x2];
x = resample(x_total, fs, fs_in);
% x = resample(sig, fs, fs_in);                                 % For demo input only

X=DFTAnaRealEntireSignal(x, K, N, p);                           % Analysis filterbank

Pxx = estimate_psd(X, lambda);
Cxx = estimate_cpsd(X(:,:,1), X(:,:,2), lambda)./ sqrt(Pxx(:,:,1).* Pxx(:,:,2));
frequency = linspace(0, fs/2, K/2+1)'; % frequency axis
Cnn = sinc(2 * frequency * d_mic/c); 

CDR = estimator(Cxx, Cnn);
CDR = max(real(CDR), 0);
CDR_avg = mean(CDR, 'all'); % average over time and frequency indices 

weights = spectral_subtraction(CDR,ss_alpha,ss_beta,mu);
weights = max(weights, G_floor);
weights = min(weights,1);
Postfilter_input = sqrt(mean(abs(X).^2,3)) .* exp(1j*angle(X(:,:,1)));
Processed = weights .* Postfilter_input;

y = DFTSynRealEntireSignal(Processed, K, N,p);                  % Synthesis filterbank


%% Plot results

figure('position',[0 0 600 250]);
imagesc(10*log10(CDR))
set(gca,'YDir','normal')
caxis([-15 15])
colorbar
title('Estimated CDR [dB]');
xlabel('Frame Index');
ylabel('Subband Index');
set(findall(gcf,'type','axes'),'fontsize',16);
set(findall(gcf,'type','text'),'fontSize',22);
% fig = gcf;
% fig.PaperPositionMode = 'auto';
% savefig(fig, 'Figures/Demo/CDR.fig');
% saveas(fig, 'Figures/Demo/CDR.png');

figure('position',[0 0 600 250]);
imagesc(weights);
set(gca,'YDir','normal');
caxis([0 1]);
colorbar;
title('Filter Gain');
xlabel('Frame Index');
ylabel('Subband Index');
set(findall(gcf,'type','axes'),'fontsize',16);
set(findall(gcf,'type','text'),'fontSize',22);
% fig = gcf;
% fig.PaperPositionMode = 'auto';
% savefig(fig, 'Figures/Demo/postfilter-gain.fig');
% saveas(fig, 'Figures/Demo/postfilter-gain.png');


