%% Description
% Evaluation of 2 Dereverberation algorithms 
% Load clean speech, RIRs and reverberant signal
% Use both algorithms to dereverb the siganl
% Evaluate both algos using the NSRR, listening test

% Reference:
% [1] Andreas Schwarz, Walter Kellermann, "Coherent-to-Diffuse Power Ratio
% Estimation for Dereverberation", IEEE/ACM Trans. on Audio, Speech and
% Lang. Proc., 2015 (under review); preprint available: arXiv:1502.03784
% PDF: http://arxiv.org/pdf/1502.03784
% [2] Nakatani T, Yoshioka T, Kinoshita K, et al. Speech Dereverberation 
% Based on Variance-Normalized Delayed Linear Prediction[J]. IEEE 
% Transactions on Audio Speech & Language Processing, 2010, 18(7):1717-1731.
% (Code by Teng Xiang at 2017-10-14)


close all;

%% Create reverberant signal

% Initial configuration
addpath(genpath('resources/Schwarz_lib'));
addpath(genpath('resources/WPE_lib'));


d_mic = 0.08;               % Mic spacing (m)
L = [6.61 5.11 2.95];       % Room dimensions [x y z] (m)
s = [1.1 2.5 1.5];          % Source position [x y z] (m) T_60 = 0.332;
r1 = [3 2.4 1.5];           % Receiver position [x y z] (m)
r2 = [3 2.4+d_mic 1.5];     % Receiver position [x y z] (m)
        
c = 342;                    % Speed of sound (m/s)
fs = 16000;                 % Sampling frequency (Hz, samples/s)
n = 4096*4;                 % Number of samples
T_60 = 0.9;                 % Reverberation time (s)

load_clean = 1;             % set to 0 if not loading clean speech siganl

if load_clean
    % Load clean speech signal s
    [sig, fs_in] = audioread('resources/clean_speech/ieee01f05.wav');

    % Create RIRs
    h1 = rir_generator(c, fs, r1, s, L, T_60, n);
    h2 = rir_generator(c, fs, r2, s, L, T_60, n);

    % Create 2-channel noisy speech
    x1 = filter(h1, 1, sig);
    x2 = filter(h2, 1, sig);
    x_total = [x1 x2];
    x = resample(x_total, fs, fs_in);

else    
    % Load an already reverberant signal (2-channel)
    % NOTE: cannot evaluate NSRR without clean ref
    [x, fs_in] =  audioread('resources/reverberant/roomC-2m-75deg.wav');
	x = resample(x, fs, fs_in);
end 



%% CDR Dereverberation

% Filterbank initialization
cdr_cfg.K = 512;                    % FFT size
cdr_cfg.N = 128;                    % Frame shift
%cdr_config.Lp = 1024;                  % Prototype filter length
load('resources/Schwarz_lib/filterbank/prototype_K512_N128_Lp1024.mat');

% Other parameters for CDR estimation
cdr_cfg.lambda = 0.68;              % Forgetting/smoothing factor for PSD estimation
cdr_cfg.mu = 1.3;                   % Noise oversubtraction/overestimation factor
cdr_cfg.G_floor = 0.1;              % Minimum gain

% Alpha and beta used for spectral subtraction
%cdr_cfg.ss_alpha = 1; cdr_cfg.ss_beta = 1;     % Power subtraction 
cdr_cfg.ss_alpha = 2; cdr_cfg.ss_beta = 0.5;    % Magnitude subtraction
%cdr_cfg.ss_alpha = 2; cdr_cfg.ss_beta = 1;     % Wiener filter

% Chosen estimator
cdr_cfg.estimator = @estimate_cdr_nodoa;   % DOA-independent estimator (CDRprop3)

% analysis filterbank
X=DFTAnaRealEntireSignal(x, cdr_cfg.K, cdr_cfg.N, p);

% estimate PSD and coherence
Pxx = estimate_psd(X, cdr_cfg.lambda);
Cxx = estimate_cpsd(X(:,:,1), X(:,:,2), cdr_cfg.lambda)./ sqrt(Pxx(:,:,1).* Pxx(:,:,2));

frequency = linspace(0, fs/2, cdr_cfg.K/2+1)'; % frequency axis

% define coherence models
Cnn = sinc(2 * frequency * d_mic/c); % diffuse noise coherence; not required for estimate_cdr_nodiffuse

% apply CDR estimator
CDR = cdr_cfg.estimator(Cxx, Cnn);
CDR = max(real(CDR), 0);
CDR_avg = mean(CDR, 'all'); % average over time and frequency indices 
CDR_avg_db = 10*log10(CDR_avg);

weights = spectral_subtraction(CDR, cdr_cfg.ss_alpha, cdr_cfg.ss_beta, cdr_cfg.mu);
weights = max(weights, cdr_cfg.G_floor);
weights = min(weights,1);

% postfilter input is computed from averaged PSDs of both microphones
Postfilter_input = sqrt(mean(abs(X).^2,3)) .* exp(1j*angle(X(:,:,1)));

% apply postfilter
Processed = weights .* Postfilter_input;

% synthesis filterbank
y_cdr = DFTSynRealEntireSignal(Processed, cdr_cfg.K, cdr_cfg.N, p);


%% WPE Dereverberation

wpe_cfg.num_mic = 2;
wpe_cfg.num_out = 2;
wpe_cfg.K = 512;                               % the number of subbands
wpe_cfg.F = 2;                                 % over-sampling rate
wpe_cfg.N = wpe_cfg.K / wpe_cfg.F;             % decimation factor
wpe_cfg.D1 = 2;                                % subband preditction delay                     
wpe_cfg.Lc1 = 30;                              % subband prediction order 
wpe_cfg.eps = 1e-4;                            % lower bound of rho(Normalizaton factor)
wpe_cfg.max_iterations = 2;

y_wpe = fdndlp(x, wpe_cfg);

%% Evaluate using NSRR

if load_clean
    ref = resample(sig, fs, fs_in);             % resample clean speech as ref
    nsrr_x = nsrr(ref, x, fs);                  % NSRR of the reverberant signal
    nsrr_y_cdr = nsrr(ref, y_cdr, fs);          % NSRR of the dereverberated siganl using CDR 
    nsrr_y_wpe = nsrr(ref, y_wpe, fs);          % NSRR of the dereverberated signal using WPE
end
%% Save audio and plot diagrams

% Save dereverberated signals
audiowrite('wav_out/reverberant.wav', x, fs);
audiowrite('wav_out/dereverberated_cdr.wav', y_cdr, fs);
audiowrite('wav_out/dereverberated_wpe.wav', y_wpe, fs);

% Plot time-frequency spectrogram
Fs = 8000;
noverlap = 128 * fs / Fs;
nfft= 256 * fs / Fs;

figure('position',[0 0 600 450]);
spectrogram(x(:,1)/max(abs(x(:,1))), hamming(nfft),noverlap,nfft,fs,'yaxis')
title('Spectrogram of reverberant signal');
set(findall(gcf,'type','axes'),'fontsize',16);
set(findall(gcf,'type','text'),'fontSize',22);
% set(gcf, 'position', [1, 235, 1366, 400]);
% set(gca, 'position', [0.05, 0.12, 0.85, 0.8]);
fig = gcf;
fig.PaperPositionMode = 'auto';
% savefig(fig, 'spectrogram-reverberant.fig');
% saveas(fig, 'spectrogram-reverberant.png');

figure('position',[0 0 600 450]);
spectrogram(y_cdr/max(abs(y_cdr)), hamming(nfft),noverlap,nfft,fs,'yaxis')
title('Spectrogram of reverberant signal');
set(findall(gcf,'type','axes'),'fontsize',16);
set(findall(gcf,'type','text'),'fontSize',22);


figure('position',[0 0 600 450]);
spectrogram(y_wpe(:,1)/max(abs(y_wpe(:,1))), hamming(nfft),noverlap,nfft,fs,'yaxis')
title('Spectrogram of reverberant signal');
set(findall(gcf,'type','axes'),'fontsize',16);
set(findall(gcf,'type','text'),'fontSize',22);

