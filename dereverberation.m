%% Description
% Evaluation of 2 Dereverberation algorithms 
% Load clean speech, RIRs and reverberant signal
% Use both algorithms to dereverb the siganl
% Evaluate both using LSD, BSD, NSRR, PESQ

% References:
% [1] Andreas Schwarz, Walter Kellermann, "Coherent-to-Diffuse Power Ratio Estimation for Dereverberation", 
% IEEE/ACM Trans. on Audio, Speech and Lang. Proc., 2015 (under review); preprint available: arXiv:1502.03784
% PDF: http://arxiv.org/pdf/1502.03784
% [2] Nakatani T, Yoshioka T, Kinoshita K, et al. Speech Dereverberation Based on Variance-Normalized Delayed 
% Linear Prediction[J]. IEEE Transactions on Audio Speech & Language Processing, 2010, 18(7):1717-1731.
% (Code by Teng Xiang at 2017-10-14)
% [3] P. A. Naylor and N. D. Gaubitch, Eds., Speech Dereverberation. Springer, 2010.
% [4] Emanuël Habets, https://www.audiolabs-erlangen.de/fau/professor/habets/software/rir-generator
% [5] ludlows, https://github.com/ludlows/pesq-mex


% Dependencies:
% 1. RIR generator
% 2. Functions from the Schwarz demo code
% 3. Functions from the WPE demo code
% 4. Evaluation measures scripts


close all;
addpath(genpath('resources/Schwarz_lib'));
addpath(genpath('resources/WPE_lib'));
addpath(genpath('resources/helper_scripts'));
addpath(genpath('resources/PESQ_lib'));

%% Create reverberant signal

% Initial configuration

% Room 1
room = 1;
L = [5.5 4 2.95];                                       % Room dimensions [x y z] (m)
T_60 = 1.5;                                             % Reverberation time (s)
% Room 2
% room = 2;
% L = [11 7 3];      
% T_60 = 3;                

d_mic = 0.08;                                           % Mic spacing (m)
c = 342;                                                % Speed of sound (m/s)
fs = 16000;                                             % Sampling frequency (Hz, samples/s)
n = 16384;                                              % Number of samples
vol = prod(L);                                          % Volume of the room (m^3)
Dc = 0.1 * sqrt( (1*vol) / (pi * T_60));                % Critical distance (m)

% FOR ROOM 1
% Config 1
config = 1;
s = [3 2 1.5];
r1 = [2.65 2 1.5];
r2 = [2.65 2+d_mic 1.5]; 

% Config 2
% config = 2;
% s = [3 2 1.5];
% r1 = [2 2 1.5];
% r2 = [2 2+d_mic 1.5]; 

% Config 3
% config = 3;
% s = [2.5 2 1.5];
% r1 = [4.5 3 1.5];
% r2 = [4.5 3+d_mic 1.5]; 

% Config 4
% config = 4;
% s = [1.5 2 1.5];
% r1 = [3 2 1.5];
% r2 = [3 2+d_mic 1.5]; 


% FOR ROOM 2
% Config 1
% config = 1;
% s = [5 3.5 1.5];
% r1 = [5.4 3.5 1.5];
% r2 = [5.4 3.5+d_mic 1.5]; 

% Config 2
% config = 2;
% s = [5 3.5 1.5];
% r1 = [7 3.5 1.5];
% r2 = [7 3.5+d_mic 1.5]; 

% Config 3
% config = 3;
% s = [5 3.5 1.5];
% r1 = [10 6 1.5];
% r2 = [10 6+d_mic 1.5]; 

% Config 4
% config = 4;
% s = [9.5 3.5 1.5];
% r1 = [5 3.5 1.5];
% r2 = [5 3.5+d_mic 1.5]; 



load_clean = 1;                                         % Set to 0 if NOT loading clean speech siganl

if load_clean
    % Load clean speech signal s
    [sig, fs_in] = audioread('resources/clean_speech/ieee01f05.wav');

    % Create RIRs
    h1 = rir_generator(c, fs, r1, s, L, T_60, n);
    h2 = rir_generator(c, fs, r2, s, L, T_60, n);

    % Create 2-channel noisy speech
    x1 = filter(h1, 1, sig);
    x1 = [x1;x1];
    x2 = filter(h2, 1, sig);
    x2 = [x2;x2];
    x_total = [x1 x2];
    x = resample(x_total, fs, fs_in);

else    
    % Load an already reverberant signal (2-channel)
    % NOTE: cannot evaluate without clean ref
    [x, fs_in] =  audioread('resources/reverberant/roomC-2m-75deg.wav');
	x = resample(x, fs, fs_in);
end 



%% CDR Dereverberation

% Filterbank initialization
cdr_cfg.K = 512;                                                        % FFT size
cdr_cfg.N = 128;                                                        % Frame shift
cdr_config.Lp = 1024;                                                   % Prototype filter length
load('resources/Schwarz_lib/filterbank/prototype_K512_N128_Lp1024.mat');

% Other parameters for CDR estimation
cdr_cfg.lambda = 0.68;                                                  % Forgetting factor for PSD estimation
cdr_cfg.mu = 1.3;                                                       % Noise oversubtraction factor
cdr_cfg.G_floor = 0.1;                                                  % Gain floor

% Alpha and beta used for spectral subtraction
%cdr_cfg.ss_alpha = 1; cdr_cfg.ss_beta = 1;                             % Power subtraction 
cdr_cfg.ss_alpha = 2; cdr_cfg.ss_beta = 0.5;                            % Magnitude subtraction
%cdr_cfg.ss_alpha = 2; cdr_cfg.ss_beta = 1;                             % Wiener filter

% Chosen estimator
cdr_cfg.estimator = @estimate_cdr_nodoa;                                % DOA-independent estimator (CDRschwarz3)

X=DFTAnaRealEntireSignal(x, cdr_cfg.K, cdr_cfg.N, p);                   % Analysis filterbank

Pxx = estimate_psd(X, cdr_cfg.lambda);
Cxx = estimate_cpsd(X(:,:,1), X(:,:,2), cdr_cfg.lambda)./ sqrt(Pxx(:,:,1).* Pxx(:,:,2));
frequency = linspace(0, fs/2, cdr_cfg.K/2+1)'; 
Cnn = sinc(2 * frequency * d_mic/c); 

% apply CDR estimator
CDR = cdr_cfg.estimator(Cxx, Cnn);
CDR = max(real(CDR), 0);
CDR_avg = mean(CDR, 'all');                                             % Averaged over time and frequency indices 
CDR_avg_db = 10*log10(CDR_avg);

weights = spectral_subtraction(CDR, cdr_cfg.ss_alpha, cdr_cfg.ss_beta, cdr_cfg.mu);
weights = max(weights, cdr_cfg.G_floor);
weights = min(weights,1);
Postfilter_input = sqrt(mean(abs(X).^2,3)) .* exp(1j*angle(X(:,:,1)));
Processed = weights .* Postfilter_input;

y_cdr = DFTSynRealEntireSignal(Processed, cdr_cfg.K, cdr_cfg.N, p);     % Synthesis filterbank


%% WPE Dereverberation

wpe_cfg.num_mic = 2;
wpe_cfg.num_out = 1;
wpe_cfg.K = 512;                                                        % Number of subbands
wpe_cfg.F = 2;                                                          % Over-sampling rate
wpe_cfg.N = wpe_cfg.K / wpe_cfg.F;                                      % Decimation factor
wpe_cfg.D1 = 2;                                                         % Subband preditction delay                     
wpe_cfg.Lc1 = 30;                                                       % Subband prediction order 
wpe_cfg.eps = 1e-4;                                                     % Lower bound of rho(Normalizaton factor)
wpe_cfg.max_iterations = 2;

y_wpe = fdndlp(x, wpe_cfg);

%% Evaluate 
if load_clean
    sig = [sig;sig];
    ref = resample(sig, fs, fs_in);             % Resample clean speech as ref

    % NSRR
    nsrr_x1 = nsrr(ref, x(:,1), fs);            % NSRR of the reverberant signal
    nsrr_x2 = nsrr(ref, x(:,2), fs);
    nsrr_y_cdr = nsrr(ref, y_cdr, fs);          % NSRR of the dereverberated siganl using CDR 
    nsrr_y_wpe = nsrr(ref, y_wpe, fs);          % NSRR of the dereverberated signal using WPE

    % Normalization
    normalized_ref = ref ./ sqrt(sum(ref .^2));
    normalized_x = x ./ sqrt(sum(x .^ 2));
    normalized_y_cdr = y_cdr ./ sqrt(sum(y_cdr .^2));
    normalized_y_wpe = y_wpe ./ sqrt(sum(y_wpe .^2));

    % LSD
    lsd_x1 = lsd(normalized_ref, normalized_x(:,1), fs);
    lsd_x2 = lsd(normalized_ref, normalized_x(:,2), fs);
    lsd_y_cdr = lsd(normalized_ref, normalized_y_cdr, fs);
    lsd_y_wpe = lsd(normalized_ref, normalized_y_wpe, fs);

    % BSD
    bsd_x1 = bsd(normalized_ref, normalized_x(:,1), fs);
    bsd_x2 = bsd(normalized_ref, normalized_x(:,2), fs);
    bsd_y_cdr = bsd(normalized_ref, normalized_y_cdr, fs);
    bsd_y_wpe = bsd(normalized_ref, normalized_y_wpe, fs);

    % PESQ
    pesq_ref = pesq_mex(normalized_ref, normalized_ref, fs, 'both');
    pesq_x1 = pesq_mex(normalized_ref, normalized_x(:,1), fs, 'both');
    pesq_x2 = pesq_mex(normalized_ref, normalized_x(:,2), fs, 'both');
    pesq_x = pesq_mex(normalized_ref, normalized_x, fs, 'both');
    pesq_y_cdr = pesq_mex(normalized_ref, normalized_y_cdr, fs, 'both');
    pesq_y_wpe = pesq_mex(normalized_ref, normalized_y_wpe, fs, 'both');


else 
    normalized_x = x ./ sqrt(sum(x .^ 2));
    normalized_y_cdr = y_cdr ./ sqrt(sum(y_cdr .^2));
    normalized_y_wpe = y_wpe ./ sqrt(sum(y_wpe .^2));
    
end


%% Save audio and plot diagrams

% Save dereverberated signals
% audiowrite('wav_out/reverberant.wav', x, fs);
% audiowrite('wav_out/dereverberated_cdr.wav', y_cdr, fs);
% audiowrite('wav_out/dereverberated_wpe.wav', y_wpe, fs);


% Plot time-frequency spectrogram
Fs = 16000;
noverlap = 128 * fs / Fs;
nfft= 256 * fs / Fs;

figure('position',[0 0 740 370]);
spectrogram(normalized_ref(:,1)/max(abs(normalized_ref(:,1))), hamming(nfft),noverlap,nfft,fs,'yaxis')
title('Spectrogram of clean signal');
set(findall(gcf,'type','axes'),'fontsize',16);
set(findall(gcf,'type','text'),'fontSize',22);
% fig = gcf;
% fig.PaperPositionMode = 'auto';
% % savefig(fig, 'Figures/Dereverb_Spec/clean.fig');
% % saveas(fig, 'Figures/Dereverb_Spec/clean.png');

figure('position',[0 0 740 370]);
spectrogram(normalized_x(:,1)/max(abs(normalized_x(:,1))), hamming(nfft),noverlap,nfft,fs,'yaxis')
title('Spectrogram of reverberant signal');
set(findall(gcf,'type','axes'),'fontsize',16);
set(findall(gcf,'type','text'),'fontSize',22);
% fig = gcf;
% fig.PaperPositionMode = 'auto';
% savefig(fig, join(['Figures/Dereverb_Spec/room', int2str(room), '-config', int2str(config), '-reverberant.fig']));
% saveas(fig, join(['Figures/Dereverb_Spec/room', int2str(room), '-config', int2str(config), '-reverberant.png']));


figure('position',[0 0 740 370]);
spectrogram(normalized_y_cdr/max(abs(normalized_y_cdr)), hamming(nfft),noverlap,nfft,fs,'yaxis')
title('Spectrogram of dereverberated signal (CDR)');
set(findall(gcf,'type','axes'),'fontsize',16);
set(findall(gcf,'type','text'),'fontSize',22);
% fig = gcf;
% fig.PaperPositionMode = 'auto';
% savefig(fig, join(['Figures/Dereverb_Spec/room', int2str(room), '-config', int2str(config), '-CDR.fig']));
% saveas(fig, join(['Figures/Dereverb_Spec/room', int2str(room), '-config', int2str(config), '-CDR.png']));


figure('position',[0 0 740 370]);
spectrogram(normalized_y_wpe/max(abs(normalized_y_wpe)), hamming(nfft),noverlap,nfft,fs,'yaxis')
title('Spectrogram of dereverberated signal (WPE)');
set(findall(gcf,'type','axes'),'fontSize',16);
set(findall(gcf,'type','text'),'fontSize',22);
% fig = gcf;
% fig.PaperPositionMode = 'auto';
% savefig(fig, join(['Figures/Dereverb_Spec/room', int2str(room), '-config', int2str(config), '-WPE.fig']));
% saveas(fig, join(['Figures/Dereverb_Spec/room', int2str(room), '-config', int2str(config), '-WPE.png']));


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
% fig = gcf;
% fig.PaperPositionMode = 'auto';
% savefig(fig, 'CDRDRR-vs-T60-positions.fig');
% saveas(fig, 'CDRDRR-vs-T60-positions.png');

figure('position',[0 0 600 450]);
plot(s(:,1),s(:,2), 'k*');
hold on
plot(r1(:,1), r1(:,2), 'b^');
plot(r2(:,1), r2(:,2), 'r^');
xlim([0 L(1)])
ylim([0 L(2)])
grid on 
legend({'Source', 'Mic 1', 'Mic 2'})
title('Positions of Source and Microphones (Top View)');
set(findall(gcf,'type','axes'),'fontsize',16);
set(findall(gcf,'type','text'),'fontSize',22);
% fig = gcf;
% fig.PaperPositionMode = 'auto';
% savefig(fig, join(['Figures/Dereverb_Config/room', int2str(room), '-config', int2str(config), '.fig']));
% saveas(fig, join(['Figures/Dereverb_Config/room', int2str(room), '-config', int2str(config), '.png']));
