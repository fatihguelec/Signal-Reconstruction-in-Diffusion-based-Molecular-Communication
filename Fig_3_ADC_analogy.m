close all
clear 
clc

fs = 2*10^2;
Ts = 1/fs;
t = 0:Ts:0.1;
fc = 10^1;
s = 2 + sin(2*pi*fc.*t);

h = figure;
plot(t, s);
hold on;
stairs(t, s);
ylim([0.9 3.1]);
ylabel('Signal Amplitude');
xlabel('Time (s)');
legend('s(t)', 'x(t)');


%% Cropped pdf figure saving
% set(h,'Units','Inches');
% pos = get(h,'Position');
% set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
% print(h,'ADC','-dpdf','-r0')



