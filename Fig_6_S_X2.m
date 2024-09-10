close all
clear
clc

%% Parameters
var_s = 10^2; %Variance of s
mu_s = 10^2; %mean value of s
T = 0.06; %sampling period
D = 10^-12; %Diffusion coefficient 
a = 1*10^-6; %radius of the NM
b = 2*10^-6; %
VRV = (4/3)*pi*b^3; %Virtual reception volume
VRN = (4/3)*pi*a^3; %Volume of Receiver NM
L = 10^6; % length


%%
s = round(normrnd(mu_s,sqrt(var_s),1,L)); %number of molecules outside
x = zeros(1,1); %Received number of molecules
temp = ((b-a)/2)+a; %Average distance of molecules to the origin
F = (a/temp).*erfc((temp-a)/sqrt(4*D*T)); %cdf
for j = 1:L
        lambda = F*s(randi(L));
%         x(j,:) = poissrnd(lambda,1,L);
        x(1,j) = poissrnd(lambda,1,1);
%         x = [x poissrnd(lambda,1,L)];
end
h = figure;
h1 = histfit(s./(VRV-VRN)); %Concentration outside
hold on;
h2 = histfit(x./VRN); %Received concentration
xlabel('Molecule Concentration (molecules/m^3)');
ylabel('Observation Frequency');
set(h1(2),'color','r'); set(h2(2),'color','g'); 
l=legend([h1(2) h2(2)],'Signal in V_R','Reconstructed Signal');
set(l,'Interpreter','tex','FontSize',8)
% title('Molecule Concentration Distribution');
delete(h1(1)); delete(h2(1));

%Find the mean and variance of the fitted distributions
pd_s = fitdist(s'./(VRV-VRN),'normal'); %Signal outside the RN in VRV
mu_1 = pd_s.mu;
var_1 = (pd_s.sigma)^2;

pd_x = fitdist(x'./VRN,'normal'); %Signal outside the RN in VRV
mu_x = pd_x.mu;
var_x = (pd_x.sigma)^2;

% set(gca,'XLim',[0 max(h2(2).XData)]);
% set(gca,'YLim',[0 1400]);
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
% print(h,'Distribution_s_x2','-dpdf','-r0')%save as pdf

% %Find maximum points
% sf = h1(2).YData;
% sf_conc = h1(2).XData;
% max_s = sf_conc(find(sf == max(sf)));
% xf = h2(2).YData;
% xf_conc = h2(2).XData;
% max_x = xf_conc(find(xf == max(xf)));
% error = max_s - max_x

