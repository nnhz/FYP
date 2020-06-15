%% Description
% CDR_bias.m
% This script is for evaluating the bias of CDR estimators

close all;

%% Initial Configuration
c = 340;
d_mic = 0.08;
f = 5000;
% TDOA = d_mic * sin(theta)/c;
TDOA = 1/(5*f);

Css = exp(1j * 2 * pi * f * TDOA);      % Gamma_s
Cnn = sinc(2 * f * d_mic/c);            % Gamma_n
CDR_log = linspace(-20,20, 201);
CDR = 10 .^ (CDR_log/10);
% theoretical CDR = (Cnn - Cxx) ./ (Cxx - Css)
Cxx = (Cnn + Css * CDR) ./ (CDR +1);    % Gamma_x

%% CDR Estimation

Jeub = max(0, (Cnn - real(conj(Css) .* Cxx)) ./ (real(conj(Css) .* Cxx) - 1));

Thiergart1 = max(0, real((Cnn - Cxx) ./ (Cxx - Css)));

Thiergart2 = max(0, real((Cnn - Cxx) ./ (Cxx - exp(1j * angle(Cxx)))));

Schwarz1 = max(0, real(conj(Css) .* (Cnn - Cxx)) ./ (real(conj(Css) .* Cxx) - 1));

Schwarz2 = max(0, 1 ./ (-abs(Cnn-exp(1j * angle(Css))) ./ (Cnn .* cos(angle(Css)) - 1)) .* abs((conj(Css) .* (Cnn - Cxx)) ./ (real(conj(Css) .* Cxx) - 1)));

Schwarz3 = max(0, (-(abs(Cxx).^2 + Cnn.^2 .* real(Cxx).^2 - Cnn.^2 .* abs(Cxx).^2 - 2 .* Cnn .* real(Cxx) + Cnn.^2).^(1/2) - abs(Cxx).^2 + Cnn .* real(Cxx)) ./ (abs(Cxx).^2 - 1));

Schwarz4 = imag(Cxx)./(imag(Css) - imag(Cxx));
Schwarz4(imag(Css)./imag(Cxx)<=1) = Inf;
Schwarz4(imag(Css)./imag(Cxx)<=0) = 0;

% Convert to the log scale
Jeub_log = 10*log10(Jeub);
Thiergart1_log = 10*log10(Thiergart1);
Thiergart2_log = 10*log10(Thiergart2);
Schwarz1_log = 10*log10(Schwarz1);
Schwarz2_log = 10*log10(Schwarz2);
Schwarz3_log = 10*log10(Schwarz3);
Schwarz4_log = 10*log10(Schwarz4);

%% Visualization of Results

figure('position',[0 0 600 450]);
plot(CDR_log, Jeub_log);
hold on;
plot(CDR_log, Thiergart1_log);
plot(CDR_log, Thiergart2_log);
plot(CDR_log, Schwarz1_log);
plot(CDR_log, Schwarz2_log);
plot(CDR_log, Schwarz3_log);
plot(CDR_log, Schwarz4_log);
hold off;
xlim([-20 20]);
ylim([-20 20]);
xlabel('Theoretical CDR/dB');
ylabel('Estimated CDR/dB');
grid on;
legend({'Jeub', 'Thiergart1', 'Thiergart2', 'Schwarz1', 'Schwarz2', 'Schwarz3', 'Schwarz4'}, 'Location', 'best');
title('CDR Estimate vs. Theoretical Value');
set(findall(gcf,'Type','Axes'),'FontSize',16);
set(findall(gcf,'Type','Text'),'FontSize',22);
set(findall(gcf,'Type','Line'),'LineWidth', 2);
% fig = gcf;
% fig.PaperPositionMode = 'auto';
% savefig(fig, join(['Figures/CDR_Evaluation/', int2str(f), 'Hz.fig']));
% saveas(fig, join(['Figures/CDR_Evaluation/', int2str(f), 'Hz.png']));
