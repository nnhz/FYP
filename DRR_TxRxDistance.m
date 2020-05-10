%% Description
% DRR_TxRxDistance.m
% Average DRR vs. Source-Receiver Distance
% Fixed config for room, source
% Vary the receiver positions (using 2 mics)


close all;

addpath(genpath('resources/Schwarz_lib'));

%% Initial parameters and configuration

L = [6.61 5.11 2.95];       % Room dimensions [x y z] (m)
s = [1.1 2.5 1.5];          % Source position [x y z] (m) T_60 = 0.332;

c = 342;                    % Speed of sound (m/s)
fs = 16000;                 % Sampling frequency (Hz, samples/s)
n = 4096*3;                 % Number of samples
T_60 = 0.45;                % Reverberation time (s)

% Number of iterations in the x- and y-direction for 2 mic positions
x_itr = 6;
y_itr = 6;

% Initialize variables for visualisation of results
distances = zeros(1, x_itr * y_itr);
r1_all = zeros(x_itr * y_itr, 3);
r2_all = zeros(x_itr * y_itr, 3);
DRR_all = zeros(x_itr * y_itr, 2);


%% Load clean speech signal
[x, fs_in] = audioread('resources/IEEE_sentences/ieee01f06.wav');


%% Obtain average DRR
% Using estDRRFromRIR from the ACE Challenge Corpus

for i = 1:y_itr
    for j = 1:x_itr

        % Generate RIRs, assuming 2 mirophones
        % https://www.audiolabs-erlangen.de/fau/professor/habets/software/rir-generator

        r1 = [1.5+j*0.5 1+i*0.5 1.5];          % Receiver position [x y z] (m)
        r2 = [1.5+j*0.5 1+d_mic+i*0.5 1.5];    % Receiver position [x y z] (m)

        r1_all(j+6*(i-1),:) = r1;
        r2_all(j+6*(i-1),:) = r2;

        distances(j+6*(i-1)) = norm(s - (r1 + r2)/2);

        h1 = rir_generator(c, fs, r1, s, L, T_60, n);
        h2 = rir_generator(c, fs, r2, s, L, T_60, n);
        
        DRR_all(j+6*(i-1), :) = estDRRFromRIR([h1;h2]', fs); % h is a vector of 2 mics?


    end

end

%% Average CDR vs. Source-Receiver Distance

DRR_avg = mean(DRR_all, 2);

figure;
scatter(distances, DRR_avg);
set(gca,'FontSize', 12);
title('Average DRR vs. Source-Receiver Distance', 'Fontsize', 18);
xlabel('Distance/m', 'Fontsize', 18);
ylabel('Average CDR', 'Fontsize', 18);
   
figure;
plot3(s(:,1), s(:,2), s(:,3), 'k*', r1_all(:,1), r1_all(:,2), r1_all(:,3), 'b^', r2_all(:,1), r2_all(:,2), r2_all(:,3), 'r>');
grid on;
legend({'s', 'r1', 'r2'})
xlim([0 L(1)])
ylim([0 L(2)])
zlim([0 L(3)])
title('Positions of Source and Microphones', 'Fontsize', 18);

