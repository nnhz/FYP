close all;

c = 340;
d_mic = 0.08;
f = 1000;
%theta = ;
%TDOA = d_mic * sin(theta)/c;
TDOA = 1/(5*f);
CDR_log = -10;


%% Cnn_diff as indepedent variable

% Cnn_diff = linspace(-0.3, 0.3, 201);
% 
% Css = exp(1j * 2 * pi * f * TDOA);              % Gamma_s
% Cnn_theo = sinc(2 * f * d_mic/c);               % Gamma_n
% Cnn = Cnn_theo + Cnn_diff;
% CDR = 10 .^ (CDR_log/10);
% Cxx = (Cnn_theo + Css * CDR) ./ (CDR +1);       % Gamma_x
% 
% Jeub = max(0, (Cnn - real(exp(-1i * angle(Css)) .* Cxx)) ./ (real(exp(-1i * angle(Css)) .* Cxx) - 1));
% 
% Thiergart1 = max(0, real((Cnn - Cxx) ./ (Cxx - Css)));
% 
% Thiergart2 = max(0, real((Cnn - Cxx) ./ (Cxx - exp(1i * angle(Cxx)))));
% 
% Schwarz1 = max(0, real(exp(-1i * angle(Css)) .* Cnn - (exp(-1i * angle(Css)) .* Cxx)) ./ (real(exp(-1i * angle(Css)) .* Cxx) - 1));
% 
% Schwarz2 = max(0, 1 ./ (-abs(Cnn-exp(1i * angle(Css))) ./ (Cnn .* cos(angle(Css)) - 1)) .* abs((exp(-1i * angle(Css)) .* Cnn - (exp(-1i * angle(Css)) .* Cxx)) ./ (real(exp(-1i * angle(Css)) .* Cxx) - 1)));
% 
% Schwarz3 = max(0, (-(abs(Cxx).^2 + Cnn.^2 .* real(Cxx).^2 - Cnn.^2 .* abs(Cxx).^2 - 2 .* Cnn .* real(Cxx) + Cnn.^2).^(1/2) - abs(Cxx).^2 + Cnn.*real(Cxx))./(abs(Cxx).^2-1));
% 
% Schwarz4 = imag(Cxx)./(imag(Css) - imag(Cxx));
% Schwarz4(imag(Css)./imag(Cxx)<=1) = Inf;
% Schwarz4(imag(Css)./imag(Cxx)<=0) = 0;
% 
% Jeub_log = 10*log10(Jeub);
% Thiergart1_log = 10*log10(Thiergart1);
% Thiergart2_log = 10*log10(Thiergart2);
% Schwarz1_log = 10*log10(Schwarz1);
% Schwarz2_log = 10*log10(Schwarz2);
% Schwarz3_log = 10*log10(Schwarz3);
% Schwarz4_log = 10*log10(Schwarz4);
% 
% Jeub_log_diff = Jeub_log - CDR_log;
% Thiergart1_log_diff = Thiergart1_log - CDR_log;
% Thiergart2_log_diff = Thiergart2_log - CDR_log;
% Schwarz1_log_diff = Schwarz1_log - CDR_log;
% Schwarz2_log_diff = Schwarz2_log - CDR_log;
% Schwarz3_log_diff = Schwarz3_log - CDR_log;
% Schwarz4_log_diff = Schwarz4_log - CDR_log;
% 
% figure('position',[0 0 800 600]);
% plot(Cnn_diff, Jeub_log_diff);
% hold on;
% plot(Cnn_diff, Thiergart1_log_diff);
% plot(Cnn_diff, Thiergart2_log_diff);
% plot(Cnn_diff, Schwarz1_log_diff);
% plot(Cnn_diff, Schwarz2_log_diff);
% plot(Cnn_diff, Schwarz3_log_diff);
% plot(Cnn_diff, ones(1, length(Cnn_diff)).* Schwarz4_log_diff); % omit as it does not depends on Cnn
% hold off;
% grid on;
% xlim([-0.3 0.3]);
% ylim([-10 10]);
% % xlabel('$$\Tilde{\Gamma}_{n} - \Gamma_{n}$$');
% xlabel('Noise Coherence Estimation Error');
% % ylabel('$$\widehat{{CDR}} - {CDR}$$ (dB)', 'Interpreter', 'Latex');
% ylabel('CDR Estimation Error');
% lgd = legend({'Jeub', 'Thiergart1', 'Thiergart2', 'Schwarz1', 'Schwarz2', 'Schwarz3', 'Schwarz4'}, 'Location', 'northeast');
% title('CDR Estimation Error vs. Noise Coherence Estimation Error');
% set(findall(gcf,'Type','Axes'),'FontSize',16);
% set(findall(gcf,'Type','Text'),'FontSize',22);
% set(findall(gcf,'Type','Line'),'LineWidth', 2);
% lgd.NumColumns = 3;
% fig = gcf;
% fig.PaperPositionMode = 'auto';
% savefig(fig, join(['Figures/CDR_Evaluation/noise-coherence-CDR', int2str(CDR_log), 'dB-f', int2str(f), 'Hz.fig']));
% saveas(fig, join(['Figures/CDR_Evaluation/noise-coherence-CDR', int2str(CDR_log), 'dB-f', int2str(f), 'Hz.png']));

%% Arg(Css) difference

Css_theo = exp(1j * 2 * pi * f * TDOA);      % Gamma_s
Cnn = sinc(2 * f * d_mic/c);            % Gamma_n
CDR = 10 .^ (CDR_log/10);
% theoretical CDR = (Cnn - Cxx) ./ (Cxx - Css)
Cxx = (Cnn + Css_theo * CDR) ./ (CDR +1);    % Gamma_x

% Css_arg_diff = linspace(-1.8, 1.8, 201);
Css_arg_diff = linspace(-2, 2, 201);
Css = Css_theo .* exp(1i * Css_arg_diff);
Css2 = wrapToPi(exp(1i*(angle(Css_theo) + Css_arg_diff)));

Jeub = max(0, (Cnn - real(exp(-1i * angle(Css)) .* Cxx)) ./ (real(exp(-1i * angle(Css)) .* Cxx) - 1));

Thiergart1 = max(0, real((Cnn - Cxx) ./ (Cxx - Css)));

Thiergart2 = max(0, real((Cnn - Cxx) ./ (Cxx - exp(1i * angle(Cxx)))));

Schwarz1 = max(0, real(exp(-1i * angle(Css)) .* Cnn - (exp(-1i * angle(Css)) .* Cxx)) ./ (real(exp(-1i * angle(Css)) .* Cxx) - 1));
% Schwarz1 =  max(0, real(conj(Css) .* (Cnn - Cxx)) ./ (real(conj(Css) .* Cxx) - 1));
% Schwarz1 = max(0, real(exp(-1i * (angle(Css_theo)+Css_arg_diff)) .* Cnn - (exp(-1i * (angle(Css_theo)+Css_arg_diff)) .* Cxx)) ./ (real(exp(-1i * (angle(Css_theo)+Css_arg_diff)) .* Cxx) - 1));


Schwarz2 = max(0, 1 ./ (-abs(Cnn-exp(1i * angle(Css))) ./ (Cnn .* cos(angle(Css)) - 1)) .* abs((exp(-1i * angle(Css)) .* Cnn - (exp(-1i * angle(Css)) .* Cxx)) ./ (real(exp(-1i * angle(Css)) .* Cxx) - 1)));

Schwarz3 = max(0, (-(abs(Cxx).^2 + Cnn.^2 .* real(Cxx).^2 - Cnn.^2 .* abs(Cxx).^2 - 2 .* Cnn .* real(Cxx) + Cnn.^2).^(1/2) - abs(Cxx).^2 + Cnn.*real(Cxx))./(abs(Cxx).^2-1));

Schwarz4 = imag(Cxx)./(imag(Css) - imag(Cxx));
Schwarz4(imag(Css)./imag(Cxx)<=1) = Inf;
Schwarz4(imag(Css)./imag(Cxx)<=0) = 0;

Jeub_log = 10*log10(Jeub);
Thiergart1_log = 10*log10(Thiergart1);
Thiergart2_log = 10*log10(Thiergart2);
Schwarz1_log = 10*log10(Schwarz1);
Schwarz2_log = 10*log10(Schwarz2);
Schwarz3_log = 10*log10(Schwarz3);
Schwarz4_log = 10*log10(Schwarz4);

Jeub_log_diff = Jeub_log - CDR_log;
Thiergart1_log_diff = Thiergart1_log - CDR_log;
Thiergart2_log_diff = Thiergart2_log - CDR_log;
Schwarz1_log_diff = Schwarz1_log - CDR_log;
Schwarz2_log_diff = Schwarz2_log - CDR_log;
Schwarz3_log_diff = Schwarz3_log - CDR_log;
Schwarz4_log_diff = Schwarz4_log - CDR_log;

figure('position',[0 0 800 600]);
plot(Css_arg_diff, Jeub_log_diff);
hold on;
plot(Css_arg_diff, Thiergart1_log_diff);
% plot(Css_arg_diff, Thiergart2_log_diff); % omitted
plot(Css_arg_diff, Schwarz1_log_diff);
plot(Css_arg_diff, Schwarz2_log_diff);
plot(Css_arg_diff, ones(1, length(Css_arg_diff)).* Schwarz3_log_diff);
plot(Css_arg_diff, Schwarz4_log_diff); 
hold off;
grid on
xlim([-1.5 1.5]);
ylim([-10 10]);
xlabel('Signal Coherence Phase Estimation Error');
% ylabel('$$\widehat{{CDR}} - {CDR}$$ (dB)', 'Interpreter', 'Latex');
ylabel('CDR Estimation Error');
lgd = legend({'Jeub', 'Thiergart1', 'Schwarz1', 'Schwarz2', 'Schwarz3', 'Schwarz4'}, 'Location', 'northeast');
title('CDR Estimation Error vs. Signal Coherence Phase Estimation Error');
set(findall(gcf,'Type','Axes'),'FontSize',16);
set(findall(gcf,'Type','Text'),'FontSize',22);
set(findall(gcf,'Type','Line'),'LineWidth', 2);
lgd.NumColumns = 3;
fig = gcf;
fig.PaperPositionMode = 'auto';
savefig(fig, join(['Figures/CDR_Evaluation/signal-coherence-phase-CDR', int2str(CDR_log), 'dB-f', int2str(f), 'Hz.fig']));
saveas(fig, join(['Figures/CDR_Evaluation/noise-coherence-phase-CDR', int2str(CDR_log), 'dB-f', int2str(f), 'Hz.png']));

