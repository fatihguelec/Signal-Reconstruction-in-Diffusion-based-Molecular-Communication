close all
clear all
clc

%% Parameters
var_s = 10^2; %Variance of s
mu_s = 10^2; %mean value of s
rho_sx = 0.75; % Correlation coeff between x and s
T = 0.1; %sampling period
D = 10^-12; %Diffusion coefficient 
a = 10^-6:0.25*10^-6:3*10^-6; %radius of the NM
% a = 10^-6:0.5*10^-6:5*10^-6;
% a = 2*10^-6:5*10^-7:6*10^-6;
b =2*a;
c1 = 1.09500814703333;
c2 = 0.75651138383854;
% yy = 0.47771*ones(1,length(a));%Roots from Mathematica for b=1.5a
yy = 0.16291*ones(1,length(a));%Roots from Mathematica for b=2a

h = figure;
linS = {'-','-','-','-.',':',':',':'};
markers = {'d','x','*','o','*','x','d'};
for i=1:length(a)
    temp = ((b(i)-a(i))/2)+a(i); %Average distance of molecules to the origin
    VRV = (4/3)*pi*b(i)^3; %Reception volume
    VRN = (4/3)*pi*a(i)^3; %Volume of Receiver NM
    F(i) = (a(i)/temp).*erfc((temp-a(i))/sqrt(4*D*T));%cdf

    % MSE for concentration with VRV and VRN with normalization of T
    MSE_sx(i) = (var_s + mu_s^2)/(VRV-VRN)^2 - 2*rho_sx*(sqrt(var_s)/(VRV-VRN))*(sqrt(F(i)*(F(i)*var_s+mu_s))/VRN)-...
    (2*F(i)*mu_s^2/(VRN*(VRV-VRN)))+((F(i)^2*var_s+F(i)*mu_s+F(i)^2*mu_s^2)/(VRN^2));    
    
    %MSE with T_opt
    z(i) = (-c1+sqrt(c1^2-4*c2*log(yy(i))))/(2*c2);
    T_opt(i) = ((b(i)-a(i))^2) / (z(i)^2*16*D);
    F_opt(i) = (a(i)/temp).*erfc((temp-a(i))/sqrt(4*D*T_opt(i)));%cdf
    MSE_opt(i) = (var_s + mu_s^2)/(VRV-VRN)^2 - 2*rho_sx*(sqrt(var_s)/(VRV-VRN))*(sqrt(F_opt(i)*(F_opt(i)*var_s+mu_s))/VRN)-...
    (2*F_opt(i)*mu_s^2/(VRN*(VRV-VRN)))+((F_opt(i)^2*var_s+F_opt(i)*mu_s+F_opt(i)^2*mu_s^2)/(VRN^2));

%     NMSE_opt(i) = MSE_opt(i)/((var_s + mu_s^2)/(VRV-VRN)^2);%normalization
    
end
NMSE_sx = MSE_sx./max(MSE_sx);%normalization
NMSE_opt = MSE_opt./max(MSE_sx);%normalization
% NMSE_sx = MSE_sx/((var_s + mu_s^2)/(VRV-VRN)^2);
plot(a, MSE_sx,'--*');
hold on;
plot(a, MSE_opt,'-o');
legend('T = 0.1 s','T = T_{opt}'); 
xlabel('a (V_N radius) (m)');
ylabel('$\mathcal{E}$(Signal Distortion)','interpreter','latex');
% xlabel('a (m)');
% ylabel('NMSE');
xlim([0.99*10^-6 3.01*10^-6]);

set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
% print(h,'Opt_T1','-dpdf','-r0')%save as pdf
