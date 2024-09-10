close all
clear all
clc

%% Parameters
var_s = 10^2; %Variance of s
mu_s = 10^2; %mean value of s
rho_sx = 0.75; % Correlation coeff between x and s
T = 0:0.01:0.25; %sampling period
% T = 0:0.01:1;

D = [10^-13 10^-12 10^-11]; %Diffusion coefficient 
% D = [5*10^-13 10^-12 2*10^-12];
a = 10^-6; %radius of the NM
b = 2*10^-6; %radius of the reception volume
VRN = (4/3)*pi*a^3; %Volume of Receiver NM
VRV = (4/3)*pi*b^3; %Reception volume

%% MSE
h = figure;
linS = {'-','-','-','-.',':',':',':'};
markers = {'d','x','*','o','*','x','d'};
temp = ((b-a)/2)+a; %Average distance of molecules to the origin
for ix = 1:1:length(D)

    for j = 1:1:length(T)
        F(j) = (a/temp).*erfc((temp-a)/sqrt(4*D(ix)*T(j)));       
       
% MSE for concentration with VRV and VRN with normalization of T
        MSE_sx(j) = (var_s + mu_s^2)/(VRV-VRN)^2 - 2*rho_sx*(sqrt(var_s)/(VRV-VRN))*(sqrt(F(j)*(F(j)*var_s+mu_s))/VRN)-...
        (2*F(j)*mu_s^2/(VRN*(VRV-VRN)))+((F(j)^2*var_s+F(j)*mu_s+F(j)^2*mu_s^2)/(VRN^2));
% MSE for concentration with VRV and VRN without normalization of T
%            MSE_sx(j) = (var_s + mu_s^2)/(VRV-VRN)^2 - 2*rho_sx(ix)*(sqrt(var_s)/(VRV-VRN))*(sqrt(F(j)^2*var_s*T(j)^2+F(j)*mu_s*T(j))/VRN)-...
%            (2*F(j)*mu_s^2*T(j)/(VRN*(VRV-VRN)))+((F(j)^2*var_s*T(j)^2+F(j)*mu_s*T(j)+F(j)^2*mu_s^2*T(j)^2)/(VRN^2));       
     
    end
%     NMSE_sx = MSE_sx/((var_s + mu_s^2)/(VRV-VRN)^2);%normalization
    plot(T, MSE_sx,'LineStyle', linS{ix}, 'Marker', markers{ix});
    legendinfo{ix} = ['D = ' num2str(D(ix)) ' m^2/s'];
    hold on;
end

xlabel('T (sampling period)(s)');
ylabel('$\mathcal{E}$(Signal Distortion)','interpreter','latex');
% xlabel('T (örnekleme periyodu)(s)');
% ylabel('NMSE');
legend(legendinfo);
% title(['\mu_s = ' num2str(mu_s) ', \sigma_s^2 = ' num2str(var_s) ', \rho_{sx} = ' num2str(rho_sx) ', b = ',num2str(b),'\mum']);
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
% print(h,'Theo_MSE_D','-dpdf','-r0')%save as pdf



