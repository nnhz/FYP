close all;

%% Initialization
n_rooms = 7;
n_mics = 2;
DRR_all = zeros(n_rooms, n_mics);
CDR_avg_all = zeros(n_rooms, 1);

% Clean input speech signal
[x, fs_in] = audioread('resources/IEEE_sentences/ieee01f06.wav');


%% Get RIRs
% Generated or from the ACE Challenge Corpus

% [h, fh_in] =  audioread(join(['resources/ACE_Corpus/Lin8Ch/RIR_', sprintf('%d', i) ,'.wav']));
[h1, fh_in] =  audioread('resources/ACE_Corpus/Lin8Ch/RIR_1.wav');  % Office 1
[h2, ~] =  audioread('resources/ACE_Corpus/Lin8Ch/RIR_2.wav');      % Office 2
[h3, ~] =  audioread('resources/ACE_Corpus/Lin8Ch/RIR_3.wav');      % Meeting Room 1
[h4, ~] =  audioread('resources/ACE_Corpus/Lin8Ch/RIR_4.wav');      % Meeting Room 2
[h5, ~] =  audioread('resources/ACE_Corpus/Lin8Ch/RIR_5.wav');      % Lecture Room 1
[h6, ~] =  audioread('resources/ACE_Corpus/Lin8Ch/RIR_6.wav');      % Lecture Room 2
[h7, ~] =  audioread('resources/ACE_Corpus/Lin8Ch/RIR_7.wav');      % Building lobby

% TODO: Refactor code

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
% use estDRRFromRIR function from the ACE Corpus
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

y1 = [filter(h1(:,1), 1, x) filter(h1(:,2),1,x)];
y2 = [filter(h2(:,1), 1, x) filter(h2(:,2),1,x)];
y3 = [filter(h3(:,1), 1, x) filter(h3(:,2),1,x)];
y4 = [filter(h4(:,1), 1, x) filter(h4(:,2),1,x)];
y5 = [filter(h5(:,1), 1, x) filter(h5(:,2),1,x)];
y6 = [filter(h6(:,1), 1, x) filter(h6(:,2),1,x)];
y7 = [filter(h7(:,1), 1, x) filter(h7(:,2),1,x)];

[CDR1, CDR_avg_all(1, :), z1] = getCDRavg(y1, fs_in);
[CDR2, CDR_avg_all(2, :), z2] = getCDRavg(y2, fs_in);
[CDR3, CDR_avg_all(3, :), z3] = getCDRavg(y3, fs_in);
[CDR4, CDR_avg_all(4, :), z4] = getCDRavg(y4, fs_in);
[CDR5, CDR_avg_all(5, :), z5] = getCDRavg(y5, fs_in);
[CDR6, CDR_avg_all(6, :), z6] = getCDRavg(y6, fs_in);
[CDR7, CDR_avg_all(7, :), z7] = getCDRavg(y7, fs_in);


%% Visualisation

% DRR vs. T_60
figure;
scatter(T_60, DRR_avg);
set(gca,'FontSize', 12);
title('Average DRR vs. T60', 'Fontsize', 18);
xlabel('T60/s', 'Fontsize', 18);
ylabel('Average DRR', 'Fontsize', 18);

% CDR vs. T_60
figure;
scatter(T_60, 10*log10(CDR_avg_all));
set(gca,'FontSize', 12);
title('Average CDR vs. T60', 'Fontsize', 18);
xlabel('T60/s', 'Fontsize', 18);
ylabel('Average CDR', 'Fontsize', 18);

% CDR vs. DRR
figure;
scatter(DRR_avg, 10*log10(CDR_avg_all));
set(gca,'FontSize', 12);
title('DRR vs. CDR', 'Fontsize', 18);
xlabel('DRR', 'Fontsize', 18);
ylabel('CDR', 'Fontsize', 18);

