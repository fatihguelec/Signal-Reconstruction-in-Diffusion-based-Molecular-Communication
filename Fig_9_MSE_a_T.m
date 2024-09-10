close all
clear
clc

%% Parameters
var_s = 10^2; %Variance of s
mu_s = 10^2; %mean value of s
rho_sx = 0.75; % Correlation coeff between x and s
T = 0:0.01:1; %sampling period
D = 10^-12; %Diffusion coefficient 
a = 1.5*10^-6:10^-7:2*10^-6; %radius of the NM
% b = a + 0.5*10^-6; %radius of the reception volume
b = 2.*a;


%% MSE
h = figure;
for k = 1:length(a)
    VRN = (4/3)*pi*a(k)^3; %Volume of Receiver NM
    temp = ((b(k)-a(k))/2)+a(k); %Average distance of molecules to the origin
    VRV = (4/3)*pi*b(k)^3; %Virtual Reception volume
    for j = 1:1:length(T)

        F(j) = (a(k)/temp).*erfc((temp-a(k))/sqrt(4*D*T(j)));%cdf
        % MSE for concentration with VRV and VRN with normalization of T
        MSE_sx(k, j) = (var_s + mu_s^2)/(VRV-VRN)^2 - 2*rho_sx*(sqrt(var_s)/(VRV-VRN))*(sqrt(F(j)*(F(j)*var_s+mu_s))/VRN)-...
        (2*F(j)*mu_s^2/(VRN*(VRV-VRN)))+((F(j)^2*var_s+F(j)*mu_s+F(j)^2*mu_s^2)/(VRN^2));
        
    end
end
NMSE_sx = MSE_sx./max(max(MSE_sx));%normalization
[MX,MY] = meshgrid(T,a);
surf(MX,MY,MSE_sx);
xlabel('T (Sampling Period) (s)');
ylabel('a (RN Radius) (m)');
zlabel('$\mathcal{E}$(Signal Distortion)','interpreter','latex');
colorbar;

% set(h,'Units','Inches');
% pos = get(h,'Position');
% set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
% print(h,'MSE_a_T','-dpdf','-r0')%save as pdf


