close all
clear
clc

%% Parameters
var_s = 10^2; %Variance of s
mu_s = 10^2; %mean value of s
D = 10^-12; %Diffusion coefficient 
a = 10^-6; %radius of the NM
% b = [2*10^-6 1.75*10^-6 1.5*10^-6];
% b = [15*10^-6 12.5*10^-6 10*10^-6]; %radius of the reception volume
b = [3*10^-6 2.5*10^-6 2*10^-6];
tau = 10^-3; %step time
delta = sqrt(D*2*tau); %step length in meters
N = 10^4; %Monte Carlo trials- You can decrease this for faster but less accurate results
VRN = (4/3)*pi*a^3; %Receiver NM volume
tp = 0:0.01:0.25;  %plot time

%% Random Walk
h = figure;
linS = {'-',':','-',':'};
markers = {'d','x','*','o'};
for ib = 1:length(b)
    VRV = (4/3)*pi*b(ib)^3; %Reception volume
    set = [1 -1];
    temp = 2;
    MSE(ib,1) = (var_s + mu_s^2)/(VRV-VRN)^2;
    for t = tp(2:length(tp)) %Simulation time in seconds
        n = t/tau; %Number of steps
        X = zeros(1, N); %number of counted molecules by NM
        sq_error = zeros(1, N);
        S = zeros(1, N); %number of remaining molecules in VRV
        for i = 1:N %Monte Carlo Loop
            for k = 1:mu_s %Loop of Molecules
                %Random initial spherical coordinates for each molecule
%                 r = (b(ib)-a).*rand(1,1) + a;
                r = ((b(ib)-a)/2)+a; %fixed initial coodinate
                theta = (pi-0).*rand(1,1);
                phi = (2*pi-0).*rand(1,1);
                %Conversion to cartesian coordinates
                x = r*sin(theta)*cos(phi);
                y = r*sin(theta)*sin(phi);
                z = r*cos(theta);

                for j = 1:n %Movement of each molecule
                    if(r <= a )
                        X(i) = X(i) + 1;
                        break;
                    end
                    x = x + delta.*set(randi(length(set)));
                    y = y + delta.*set(randi(length(set)));
                    z = z + delta.*set(randi(length(set)));
                    r = sqrt(x^2+y^2+z^2);
                end
            end
            S(i) = mu_s-X(i);
            sq_error(i) = ((S(i)/(VRV-VRN)) - (X(i)/VRN))^2 ;     
        end
        X_mean(temp) = mean(X);
        MSE(ib, temp) = mean(sq_error);
        temp = temp + 1;
    end
    plot(tp, MSE(ib,:),'LineStyle', linS{ib}, 'Marker', markers{ib}); %unnormalized    
    legendinfo{ib} = ['b = ' num2str(b(ib)) 'm'];
    hold on;
end

xlabel('Time(s)');
ylabel('$\mathcal{E}$(Signal Distortion)','interpreter','latex');
legend(legendinfo);
title('Random Walk Simulation');
