%% Description
% DRR_TxRxDistance.m
% Average DRR vs. Source-Receiver Distance
% Fixed config for room, source, receiver separation
% Vary the microphone (receiver) positions (using 2 mics)
% Average DRR obtained across all frequency and time indices

% References:
% [1] P. A. Naylor and N. D. Gaubitch, Eds., Speech Dereverberation. Springer, 2010.
% [2] Emanuël Habets, https://www.audiolabs-erlangen.de/fau/professor/habets/software/rir-generator

% Dependencies:
% 1. RIR generator
% 2. estDRRFromRIR (from [1])

close all;

addpath(genpath('resources/helper_scripts'));

%% Initial parameters and configuration

L = [6.61 5.11 2.95];                       % Room dimensions [x y z] (m)
s = [1.1 2.5 1.5];                          % Source position [x y z] (m) T_60 = 0.332;

c = 342;                                    % Speed of sound (m/s)
fs = 16000;                                 % Sampling frequency (Hz, samples/s)
n = 12288;                                  % Number of samples
T_60 = 0.6;                                 % Reverberation time (s)
d_mic = 0.06;                               % Mic spacing (m)

% Number of iterations in the x- and y-direction for 2 mic positions
x_itr = 6;
y_itr = 6;

% Initialize variables for visualisation of results
distances = zeros(1, x_itr * y_itr);
r1_all = zeros(x_itr * y_itr, 3);
r2_all = zeros(x_itr * y_itr, 3);
DRR_all = zeros(x_itr * y_itr, 2);

%% Obtain average DRR

for i = 1:y_itr
    for j = 1:x_itr

        r1 = [1.5+j*0.5 1+i*0.5 1.5];          % Mic 1 position [x y z] (m)
        r2 = [1.5+j*0.5 1+d_mic+i*0.5 1.5];    % Mic 2 position [x y z] (m)

        r1_all(j+6*(i-1),:) = r1;
        r2_all(j+6*(i-1),:) = r2;

        distances(j+6*(i-1)) = norm(s - (r1 + r2)/2);

        h1 = rir_generator(c, fs, r1, s, L, T_60, n);
        h2 = rir_generator(c, fs, r2, s, L, T_60, n);
        
        DRR_all(j+6*(i-1), :) = estDRRFromRIR([h1;h2]', fs); 
    end

end

%% Average DRR vs. Source-Receiver Distance

DRR_avg = mean(DRR_all, 2);

figure('position',[0 0 600 450]);
scatter(distances, DRR_avg, 80, 'filled');
title('Average DRR vs. Source-Receiver Distance');
xlabel('Distance/m');
ylabel('Average DRR/dB');
grid on;
set(findall(gcf,'type','axes'),'fontsize',16);
set(findall(gcf,'type','text'),'fontSize',22);
% fig = gcf;
% fig.PaperPositionMode = 'auto';
% savefig(fig, 'DRR-vs-TxRxDistance.fig');
% saveas(fig, 'DRR-vs-TxRxDistance.png');

figure('position',[0 0 600 450]);
plot3(s(:,1), s(:,2), s(:,3), 'k*', r1_all(:,1), r1_all(:,2), r1_all(:,3), 'b^', r2_all(:,1), r2_all(:,2), r2_all(:,3), 'r>');
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
% savefig(fig, 'DRR-vs-TxRxDistance-positions.fig');
% saveas(fig, 'DRR-vs-TxRxDistance-positions.png');

