close all;

c = 340;
d_mic = 0.06;
f = 5000;
%theta = ;
%TDOA = d_mic * sin(theta)/c;
TDOA = 1/(5*f);

Css = exp(1j * 2 * pi * f * TDOA);  
Cnn = sinc(2 * f * d_mic/c); 

CDR_log = linspace(-20,20, 201);
% theoretical CDR = (Cnn - Cxx) ./ (Cxx - Css)
CDR = 10 .^ (CDR_log/10);
Cxx = (Cnn + Css * CDR) ./ (CDR +1);


Jeub = max(0, (Cnn - real(exp(-1i*angle(Css)).*Cxx)) ./ (real(exp(-1i*angle(Css)).*Cxx) - 1));

Thiergart1 = max(0, real((Cnn - Cxx) ./ (Cxx - Css)));

Thiergart2 = max(0, real((Cnn - Cxx) ./ (Cxx - exp(1i*angle(Cxx)))));

Schwarz1 = max(0, real(exp(-1j*angle(Css)).*Cnn - (exp(-1i*angle(Css)).*Cxx))./(real(exp(-1i*angle(Css)).*Cxx) - 1));

Schwarz2 = max(0, 1./(-abs(Cnn-exp(1j*angle(Css)))./(Cnn.*cos(angle(Css))-1)).*abs((exp(-1j*angle(Css)).*Cnn - (exp(-1i*angle(Css)).*Cxx))./(real(exp(-1i*angle(Css)).*Cxx) - 1)));

Schwarz3 = max(0, (-(abs(Cxx).^2 + Cnn.^2.*real(Cxx).^2 - Cnn.^2.*abs(Cxx).^2 - 2.*Cnn.*real(Cxx) + Cnn.^2).^(1/2) - abs(Cxx).^2 + Cnn.*real(Cxx))./(abs(Cxx).^2-1));

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


figure;
plot(CDR_log, Jeub_log);
hold on;
plot(CDR_log, Thiergart1_log);
plot(CDR_log, Thiergart2_log);
plot(CDR_log, Schwarz1_log);
plot(CDR_log, Schwarz2_log);
plot(CDR_log, Schwarz3_log);
plot(CDR_log, Schwarz4_log);
hold off;
xlabel('Theoretical CDR/dB');
ylabel('Estimated CDR/dB');
legend({'Jeub', 'Thiergart1', 'Thiergart2', 'Schwarz1', 'Schwarz2', 'Schwarz3', 'Schwarz4'});
title(join(['f = ', int2str(f), ' Hz']));
set(findall(gcf,'type','axes'),'fontsize',16);
set(findall(gcf,'type','text'),'fontSize',22);
% fig = gcf;
% fig.PaperPositionMode = 'auto';
% savefig(fig, 'CDR-vs-TxRxDistance-positions.fig');
% saveas(fig, 'CDR-vs-TxRxDistance-positions.png');
