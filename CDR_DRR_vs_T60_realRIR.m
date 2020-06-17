%% Description
% CDR_DRR_vs_T60_realRIR.m
% To investigate the relationship between CDR and T60 using real RIRs from the ACE Corpus

% References:
% [1] Andreas Schwarz, Walter Kellermann, "Coherent-to-Diffuse Power Ratio Estimation for Dereverberation", 
% IEEE/ACM Trans. on Audio, Speech and Lang. Proc., 2015 (under review); preprint available: arXiv:1502.03784
% PDF: http://arxiv.org/pdf/1502.03784
% [2] Emanuël Habets, https://www.audiolabs-erlangen.de/fau/professor/habets/software/rir-generator
% [3] P. A. Naylor and N. D. Gaubitch, Eds., Speech Dereverberation. Springer, 2010.
% [4] J. Eaton, N. D. Gaubitch, A. H. Moore, and P. A. Naylor, "The ACE challenge - Corpus Description 
% and Performance Evaluation", in 2015 IEEE Workshop on Applications of Signal Processing to Audio and Acoustics 
% (WASPAA), pp. 1-5, Oct 2015

% Dependencies:
% 1. RIR generator
% 2. Functions from the Schwarz demo code (wrapped in getCDRavg.m)
% 3. estDRRFromRIR
% 4. Real RIRs from the ACE Corpus

close all;

addpath(genpath('resources/helper_scripts'));

%% Initialization
n_rooms = 7;
n_mics = 2;
DRR_all = zeros(n_rooms, n_mics);
CDR_avg_all = zeros(n_rooms, 1);

% Clean input speech signal
[sig, fs_in] = audioread('resources/clean_speech/ieee01f05.wav');


%% Get RIRs
% Lengths of RIRs are different
[h1, fh_in] =  audioread('resources/ACE_Corpus/Lin8Ch/RIR_1.wav');  % Office 1
[h2, ~] =  audioread('resources/ACE_Corpus/Lin8Ch/RIR_2.wav');      % Office 2
[h3, ~] =  audioread('resources/ACE_Corpus/Lin8Ch/RIR_3.wav');      % Meeting Room 1
[h4, ~] =  audioread('resources/ACE_Corpus/Lin8Ch/RIR_4.wav');      % Meeting Room 2
[h5, ~] =  audioread('resources/ACE_Corpus/Lin8Ch/RIR_5.wav');      % Lecture Room 1
[h6, ~] =  audioread('resources/ACE_Corpus/Lin8Ch/RIR_6.wav');      % Lecture Room 2
[h7, ~] =  audioread('resources/ACE_Corpus/Lin8Ch/RIR_7.wav');      % Building lobby

% Only use 2 mics (of the 8 mics)
h1 = h1(:, 4:5);
h2 = h2(:, 4:5);
h3 = h3(:, 4:5);
h4 = h4(:, 4:5);
h5 = h5(:, 4:5);
h6 = h6(:, 4:5);
h7 = h7(:, 4:5);

T_60 = [0.332; 0.39; 0.437; 0.371; 0.638; 1.22; 0.646];

%% Get DRR
% Output from estDRRFromRIR is already in dB

[ DRR_all(1,:), ~ ] = estDRRFromRIR(h1, fh_in);
[ DRR_all(2,:), ~ ] = estDRRFromRIR(h2, fh_in);
[ DRR_all(3,:), ~ ] = estDRRFromRIR(h3, fh_in);
[ DRR_all(4,:), ~ ] = estDRRFromRIR(h4, fh_in);
[ DRR_all(5,:), ~ ] = estDRRFromRIR(h5, fh_in);
[ DRR_all(6,:), ~ ] = estDRRFromRIR(h6, fh_in);
[ DRR_all(7,:), ~ ] = estDRRFromRIR(h7, fh_in);

DRR_avg = mean(DRR_all, 2);

%% Get CDR estimate
y1 = [filter(h1(:,1), 1, sig) filter(h1(:,2),1,sig)];
y2 = [filter(h2(:,1), 1, sig) filter(h2(:,2),1,sig)];
y3 = [filter(h3(:,1), 1, sig) filter(h3(:,2),1,sig)];
y4 = [filter(h4(:,1), 1, sig) filter(h4(:,2),1,sig)];
y5 = [filter(h5(:,1), 1, sig) filter(h5(:,2),1,sig)];
y6 = [filter(h6(:,1), 1, sig) filter(h6(:,2),1,sig)];
y7 = [filter(h7(:,1), 1, sig) filter(h7(:,2),1,sig)];

[CDR1, CDR_avg_all(1, :), z1] = getCDRavg(y1, fs_in);
[CDR2, CDR_avg_all(2, :), z2] = getCDRavg(y2, fs_in);
[CDR3, CDR_avg_all(3, :), z3] = getCDRavg(y3, fs_in);
[CDR4, CDR_avg_all(4, :), z4] = getCDRavg(y4, fs_in);
[CDR5, CDR_avg_all(5, :), z5] = getCDRavg(y5, fs_in);
[CDR6, CDR_avg_all(6, :), z6] = getCDRavg(y6, fs_in);
[CDR7, CDR_avg_all(7, :), z7] = getCDRavg(y7, fs_in);

%% Visualisation

% DRR vs. T_60
figure('position',[0 0 600 450]);
scatter(T_60, DRR_avg, 80, 'filled');
title('Average DRR vs. T60');
xlabel('T60/s');
ylabel('Average DRR/dB');
grid on;
set(findall(gcf,'type','axes'),'fontsize',16);
set(findall(gcf,'type','text'),'fontSize',22);
% fig = gcf;
% fig.PaperPositionMode = 'auto';
% savefig(fig, 'DRR-vs-T60.fig');
% saveas(fig, 'DRR-vs-T60.png');


% CDR vs. T_60
figure('position',[0 0 600 450]);
scatter(T_60, 10*log10(CDR_avg_all), 80, 'filled');
title('Average CDR vs. T60');
xlabel('T60/s');
ylabel('Average CDR/dB');
grid on;
set(findall(gcf,'type','axes'),'fontsize',16);
set(findall(gcf,'type','text'),'fontSize',22);
% fig = gcf;
% fig.PaperPositionMode = 'auto';
% savefig(fig, 'CDR-vs-T60.fig');
% saveas(fig, 'CDR-vs-T60.png');


% CDR vs. DRR
figure('position',[0 0 600 450]);
scatter(DRR_avg, 10*log10(CDR_avg_all), 80, 'filled');
title('Average CDR vs. Average DRR');
xlabel('Average DRR/dB');
ylabel('Average CDR/dB');
grid on;
set(findall(gcf,'type','axes'),'fontsize',16);
set(findall(gcf,'type','text'),'fontSize',22);
% fig = gcf;
% fig.PaperPositionMode = 'auto';
% savefig(fig, 'CDR-vs-DRR.fig');
% saveas(fig, 'CDR-vs-DRR.png');

