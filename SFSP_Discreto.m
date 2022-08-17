clear; close all; clc

s = tf('s');
Ts = 0.5;
z = tf('z',Ts);

%config disturbance parameters (sine)
A_sin=0.6; Freq=1.35/(2*pi); 

%config disturbance parameters (step)
A_step=1e-2; t_step=50;

% choose s=1 for sine disturbance or s=0 for step disturbance
sd=0;


% Caso nominal
Tsim = 100;
Truido = 100;
ref=1;
NoisePower = 1e-7;


%% Modelos contínuos MIMO
G = 1/(s^2);
L = 5;

% Modelo do processo
Gr = G;
Lr = L; % caso nominal

%% Modelo nominal discreto
d = round(L/Ts);
Gd = c2d(G,Ts);
[nG,dG] =  tfdata(Gd,'v');

%polos para o sistema de equações para V
p_i = pole(Gd)';

%% filtro de robustez
wdist=0.1;
ppp=roots(-[1 -2*cos(wdist*Ts) 1])'; % perturbação seno


bV=[0 0 0];
pV=[p_i 1];
p=[0 0];
% return
% [Fr, kr, V, V_, S_,S1,w]=ProjSFSP_d(A,B,C,D,Ts,Lr,p,Beta,P_Rej);
[Fr, kr, V, V_, S_]=SFSPd(Gd,Ts,Lr,p,bV,pV);
S=minreal(S_-V_*z^(-d),1e-2);


%% Simulations
sim('sim_sfsp_Discreto');
figure(1)
subplot(2,1,1); %% separa a figura 1 em 2 linhas e escolhe a 1°
%plotagens
plot(t,r,'--','LineWidth',1.5,'Color','Green')%referencia
hold on
plot(t,y,'LineWidth',1,'Color','Blue')%referencia
axis([min(t) max(t) 0 max(y)*1.2])

subplot(2,1,2); %% escolhe a 2°
plot(t,u,'LineWidth',1,'Color','Yellow')%referencia
axis([min(t) max(t) min(u) max(u)*1.2])










