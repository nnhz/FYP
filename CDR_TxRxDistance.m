%% Description
% CDR_TxRxDistance.m
% Average CDR vs. Source-Receiver Distance
% Fixed config for room, source, receiver separation
% Vary the microphone (receiver) positions (using 2 mics)
% Average CDR obtained across all frequency and time indices

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
vol = prod(L);                                              % Volume of the room (m^3)
c = 342;                                                    % Speed of sound (m/s)
fs = 16000;                                                 % Sampling frequency (Hz, samples/s)
n = 12288;                                                  % Number of samples
T_60 = 0.6;                                                 % Reverberation time (s)
d_mic = 0.06;                                               % Mic spacing (m)
Dc= 0.1 * sqrt( (1*vol) / (pi * T_60));                     % Critical distance (m)

% Number of iterations in the x- and y-direction for 2 mic positions
x_itr = 6;
y_itr = 6;

% Initialize variables for results visualization
CDR_avg = zeros(1, x_itr * y_itr);
distances = zeros(1, x_itr * y_itr);
r1_all = zeros(x_itr * y_itr, 3);
r2_all = zeros(x_itr * y_itr, 3);

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

for i = 1:y_itr
    for j = 1:x_itr

        % Generate RIRs, assuming 2 mirophones
        % https://www.audiolabs-erlangen.de/fau/professor/habets/software/rir-generator

        r1 = [1.5+j*0.5 1+i*0.5 1.5];                       % Mic 1 position [x y z] (m)
        r2 = [1.5+j*0.5 1+d_mic+i*0.5 1.5];                 % Mic 2 position [x y z] (m)

        r1_all(j+6*(i-1),:) = r1;
        r2_all(j+6*(i-1),:) = r2;

        distances(j+6*(i-1)) = norm(s - (r1 + r2)/2);

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
        CDR_avg(j+6*(i-1)) = mean(CDR, 'all');              % Averaged over time and frequency indices 
        
        weights = spectral_subtraction(CDR,ss_alpha,ss_beta,mu);
        weights = max(weights, G_floor);
        weights = min(weights,1);
        Postfilter_input = sqrt(mean(abs(X).^2,3)) .* exp(1j*angle(X(:,:,1)));
        Processed = weights .* Postfilter_input;

        y = DFTSynRealEntireSignal(Processed, K, N,p);      % Synthesis filterbank

    end

end

%% Average CDR vs. Source-Receiver Distance
   
figure('position',[0 0 600 450]);
scatter(distances, 10*log10(CDR_avg), 80, 'filled');
title('Average CDR vs. Source-Receiver Distance');
xlabel('Distance/m');
ylabel('Average CDR/dB');
grid on;
set(findall(gcf,'type','axes'),'fontsize',16);
set(findall(gcf,'type','text'),'fontSize',22);
% fig = gcf;
% fig.PaperPositionMode = 'auto';
% savefig(fig, 'CDR-vs-TxRxDistance.fig');
% saveas(fig, 'CDR-vs-TxRxDistance.png');

figure('position',[0 0 600 450]);
plot3(s(:,1), s(:,2), s(:,3), 'k*', r1_all(:,1), r1_all(:,2), r1_all(:,3), 'b^', r2_all(:,1), r2_all(:,2), r2_all(:,3), 'r>');
grid on;
legend({'Source', 'Mic 1', 'Mic 2'})
xlim([0 L(1)])
ylim([0 L(2)])
zlim([0 L(3)])
xlabel('x');
ylabel('y');
zlabel('z');
title('Positions of Source and Microphones');
set(findall(gcf,'type','axes'),'fontsize',16);
set(findall(gcf,'type','text'),'fontSize',22);
% fig = gcf;
% fig.PaperPositionMode = 'auto';
% savefig(fig, 'CDR-vs-TxRxDistance-positions.fig');
% saveas(fig, 'CDR-vs-TxRxDistance-positions.png');


