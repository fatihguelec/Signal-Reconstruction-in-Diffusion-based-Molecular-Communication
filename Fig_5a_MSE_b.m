close all
clear
clc

%% Parameters
var_s = 10^2; %Variance of s
mu_s = 10^2; %mean value of s
rho_sx = 0.75; % Correlation coeff between x and s
T = 0:0.01:0.25; %sampling period
D = 10^-12; %Diffusion coefficient 
a = 10^-6; %radius of the NM
% b = [2*10^-6 1.75*10^-6 1.5*10^-6]; %radius of the reception volume
b = [3*10^-6 2.5*10^-6 2*10^-6];
% b = [15*10^-6 12.5*10^-6 10*10^-6]; %radius of the reception volume
VRN = (4/3)*pi*a^3; %Volume of Receiver NM

%% MSE
h = figure;
linS = {'-','-','-','-.',':',':',':'};
markers = {'d','x','*','o','*','x','d'};

for ix = 1:1:length(b)
    temp = ((b(ix)-a)/2)+a; %Average distance of molecules to the origin
    VRV = (4/3)*pi*b(ix)^3; %Reception volume
    for j = 1:1:length(T)

        F(j) = (a/temp).*erfc((temp-a)/sqrt(4*D*T(j)));%cdf

% MSE for concentration with VRV and VRN with normalization of T
           MSE_sx(ix, j) = (var_s + mu_s^2)/(VRV-VRN)^2 - 2*rho_sx*(sqrt(var_s)/(VRV-VRN))*(sqrt(F(j)*(F(j)*var_s+mu_s))/VRN)-...
           (2*F(j)*mu_s^2/(VRN*(VRV-VRN)))+((F(j)^2*var_s+F(j)*mu_s+F(j)^2*mu_s^2)/(VRN^2));
    end

end

for ix = 1:1:length(b)
    plot(T, MSE_sx(ix,:),'LineStyle', linS{ix}, 'Marker', markers{ix});
    legendinfo{ix} = ['b = ' num2str(b(ix)) ' m'];
    hold on;
end

xlabel('T (sampling period)(s)');
ylabel('$\mathcal{E}$(Signal Distortion)','interpreter','latex');
legend(legendinfo);
title('Theoretical Calculation');

% set(h,'Units','Inches');
% pos = get(h,'Position');
% set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
% print(h,'MSE_b3','-dpdf','-r0')%save as pdf




