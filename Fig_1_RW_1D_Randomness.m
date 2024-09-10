close all
clear all
clc

%% Parameters
Q = 10^4; %number of molecules
a = 10^-6;%radius of RN (m)
d = 10*10^-6; %distance of molecules to TN (m)
t = 0:0.05:10; %sampling period (s)
D = 10^-11; %Diffusion coefficient (m^2/s)

h = figure;
%%Random Walk 1-D
tau = 10^-3; %step time(s)
delta = sqrt(D*2*tau); %step length in m
setr = [1 -1];
temp = 1;
S = zeros(1, length(t));%initial number of molecules in V
for i = 1:length(t) %Simulation time in seconds
    n = t(i)/tau; %Number of steps
        for k = 1:Q %Loop of Molecules
            %initial coordinate for each molecule
            x = 0;
            for j = 1:n %Movement of each molecule
                x = x + delta.*setr(randi(length(setr)));            
            end
            if abs(d-x) <= a
                S(i) = S(i) + 1;
            end
        end
    C_rw(i) = S(i)/(2*a);  %concentration 
end
plot(t, C_rw);
hold on;

%Theoretical
for i = 1:length(t)
    C(i) = (Q/sqrt(4*pi*D*t(i)))*exp(-d^2/(4*D*t(i))); %single puff emission
end

plot(t, C, 'LineWidth',2);
xlabel('Time (s)');
ylabel('Molecule Concentration (molecules/m)');
legend('1-D Random Walk Simulation','Theoretical Model in (1)');

set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
print(h,'Proof_Random','-dpdf','-r0')
