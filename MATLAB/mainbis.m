clear all; close all; clc;
%%
% Constants
Re = 7.9; % ohm
Le0 = 5.65*10^(-3); % H
Bl0 = 14.54; % T.m
k0 = 6259.3; % N/m
Rm = 2.628; % N.s/m
mt = 52.08*10^(-3); % kg
b1 = 2.6826; % T
b2 = -4.0269*10^3; % T/m
k1 = -1.4562*10^5; % N/m2
k2 = 2.6106*10^7; % N/m3
l1 = -5.2937*10^(-1); % H/m
l2 = -1.2012*10^2; % H/m2

% Simulation
TIME_SIM = 5;
STEP_SIZE = 0.0001;
t=0:STEP_SIZE:TIME_SIM;

% P3
Au = 5; % V
fc = 20; % Hz
ue = Au*sin(2*pi*fc*t);

[tsp3,xsp3,ys] = sim('modelP3',TIME_SIM,[],[t' ue']);
uer = ue(1/STEP_SIZE:length(t));  %uer = ue; % uer = ue(1/STEP_SIZE:length(t)); 
xsrp3 = xsp3(1/STEP_SIZE:length(t),:); %xsr = xs; % xsr = xs(1/STEP_SIZE:length(t),:);
tr = t(1/STEP_SIZE:length(t)); %tr = t; % tr = t(1/STEP_SIZE:length(t));

%%
figure
axP = axes; set(axP, 'FontSize', 14)
grid on,subplot(411), plot(t(1:2/STEP_SIZE),ue(1:2/STEP_SIZE))
title('input $u_e$ and states responses in time domain', 'FontSize', 14, 'Interpreter','Latex')
xlabel('Time [s]', 'FontSize', 14, 'Interpreter','Latex')
ylabel('$u_e$ [V]', 'FontSize', 14, 'Interpreter','Latex')
grid on,subplot(412), plot(t(1:2/STEP_SIZE),xsp3((1:2/STEP_SIZE),1))
xlabel('Time [s]', 'FontSize', 14, 'Interpreter','Latex')
ylabel('$x$ [m]', 'FontSize', 14, 'Interpreter','Latex')
grid on,subplot(413), plot(t(1:2/STEP_SIZE),xsp3((1:2/STEP_SIZE),2))
xlabel('Time [s]', 'FontSize', 14, 'Interpreter','Latex')
ylabel('$\dot{x}$ [m/s]', 'FontSize', 14, 'Interpreter','Latex')
grid on,subplot(414), plot(t(1:2/STEP_SIZE),xsp3((1:2/STEP_SIZE),3))
xlabel('Time [s]', 'FontSize', 14, 'Interpreter','Latex')
ylabel('$i$ [A]', 'FontSize', 14, 'Interpreter','Latex')

%%
Fs=1/STEP_SIZE;
Pxxp3=[];
fp3=[];
for i = 1:1:3
    [Pxx_tmp, f_tmp]=power_spectral_density(xsrp3(:,i), Fs);
    Pxxp3=[Pxxp3, Pxx_tmp];
    fp3=[fp3, f_tmp'];
end
[UE, fe]=power_spectral_density(uer, Fs);

Pxxp3_db=10*log10(Pxxp3);
UE_db=10*log10(UE);

figure
grid on, axP = axes; set(axP, 'FontSize', 14)
subplot(411), plot(fe,UE_db)
title('Periodogram Using FFT')
axis([fc-10 5*fc -inf inf]);
xlabel('frequency [Hz]', 'FontSize', 14, 'Interpreter','Latex')
ylabel('$P_{u_e}$ [dB/Hz]', 'FontSize', 14, 'Interpreter','Latex')
subplot(412), plot(fp3(:,1),Pxxp3_db(:,1))
axis([fc-10 5*fc -inf inf]);
xlabel('frequency [Hz]', 'FontSize', 14, 'Interpreter','Latex')
ylabel('$P_{x}$ [dB/Hz]', 'FontSize', 14, 'Interpreter','Latex')
subplot(413), plot(fp3(:,2),Pxxp3_db(:,2))
axis([fc-10 5*fc -inf inf]);
xlabel('frequency [Hz]', 'FontSize', 14, 'Interpreter','Latex')
ylabel('$P_{\dot{x}}$ [dB/Hz]', 'FontSize', 14, 'Interpreter','Latex')
subplot(414), plot(fp3(:,3),Pxxp3_db(:,3))
axis([fc-10 5*fc -inf inf]);
xlabel('frequency [Hz]', 'FontSize', 14, 'Interpreter','Latex')
ylabel('$P_{i}$ [dB/Hz]', 'FontSize', 14, 'Interpreter','Latex')

%% PB4

Amplitudep4=[];
for j=1:3
    for i=1:6
       Amplitudep4(j,i)=amplituder(xsrp3(:,j), TIME_SIM, fc*i);
    end
end

NumberHarmonic = 6;
THDpb4 = THD(Amplitudep4(2,:),NumberHarmonic);
d2 = Amplitudep4(2,2)/sqrt(Amplitudep4(2,1)^2+Amplitudep4(2,2)^2)*100;
d3 = Amplitudep4(2,3)/sqrt(Amplitudep4(2,1)^2+Amplitudep4(2,3)^2)*100;

THDMatlab = thd(xsrp3(:,2), Fs, NumberHarmonic);

%% PB5

% AuPb5 = 2.5:2.5:10; % V
% fcPb5 = 20:5:200; % Hz
% d2Pb5 = zeros(length(AuPb5),length(fcPb5));
% d3Pb5 = zeros(length(AuPb5),length(fcPb5));
% 
% for k=1:length(AuPb5)
%     for l=1:length(fcPb5)
%         uePb5 = AuPb5(k)*sin(2*pi*fcPb5(l)*t);
%         [tsPb5,xsPb5,ysPb5] = sim('modelP3',TIME_SIM,[],[t' uePb5']);
%         xsrPb5 = xsPb5(1/STEP_SIZE:length(t),:);
%         AmplitudePb5=[];
%         for j=1:3
%             for i=1:3
%                AmplitudePb5(j,i)=amplituder(xsrPb5(:,j), TIME_SIM, fcPb5(l)*i);
%             end
%         end
%         d2Pb5(k,l) = AmplitudePb5(2,2)/sqrt(AmplitudePb5(2,1)^2+AmplitudePb5(2,2)^2)*100;
%         d3Pb5(k,l) = AmplitudePb5(2,3)/sqrt(AmplitudePb5(2,1)^2+AmplitudePb5(2,3)^2)*100;
%     end
% end
% 
% figure
% grid on, axP = axes; set(axP, 'FontSize', 14)
% subplot(411)
% plot(fcPb5, d2Pb5(1,:))
% hold on
% plot(fcPb5, d3Pb5(1,:))
% xlabel('Frequency [Hz]', 'FontSize', 14, 'Interpreter','Latex')
% ylabel('$A_u = 2.5$ V', 'FontSize', 14, 'Interpreter','Latex')
% l=legend('$d_2 [\%]$','$d_3 [\%]$');
% set(l, 'FontSize',14, 'Interpreter','Latex')
% subplot(412)
% plot(fcPb5, d2Pb5(2,:))
% hold on
% plot(fcPb5, d3Pb5(2,:))
% xlabel('Frequency [Hz]', 'FontSize', 14, 'Interpreter','Latex')
% ylabel('$A_u = 5$ V', 'FontSize', 14, 'Interpreter','Latex')
% l=legend('$d_2 [\%]$','$d_3 [\%]$');
% set(l, 'FontSize',14, 'Interpreter','Latex')
% subplot(413)
% plot(fcPb5, d2Pb5(3,:))
% hold on
% plot(fcPb5, d3Pb5(3,:))
% xlabel('Frequency [Hz]', 'FontSize', 14, 'Interpreter','Latex')
% ylabel('$A_u = 7.5$ V', 'FontSize', 14, 'Interpreter','Latex')
% l=legend('$d_2 [\%]$','$d_3 [\%]$');
% set(l, 'FontSize',14, 'interpreter','Latex')
% subplot(414)
% plot(fcPb5, d2Pb5(4,:))
% hold on
% plot(fcPb5, d3Pb5(4,:))
% xlabel('Frequency [Hz]', 'FontSize', 14, 'Interpreter','Latex')
% ylabel('$A_u = 10$ V', 'FontSize', 14, 'Interpreter','Latex')
% l=legend('$d_2 [\%]$','$d_3 [\%]$');
% set(l, 'interpreter','Latex', 'FontSize',14)

%% P6

[x,u,y,dx,options] = trim('modelP3',[0;0;0],[0],[],[],[1],[]);

[A,B,C,D] = linmod('modelP3',[0;0;0],0);


%% P7

Mc = [B A*B A^2*B];
rank(Mc); % = 3 controllable
Mo=[C
    C*A
    C*A^2];
rank(Mo); % = 3 observable


%% P8

lambda = eig(A); % lambda in the left half plane -> stable
syms X
eq = X^3 + (Re/Le0 + Rm/mt)*X^2 + (Rm*Re+Bl0^2+k0*Le0)/(mt*Le0)*X + k0*Re/(mt*Le0);

eval(subs(eq, X, lambda(1)));
eval(subs(eq, X, lambda(2)));
eval(subs(eq, X, lambda(3)));

%% PB9
%ue = Au*sin(2*pi*fc*t) + Au*sin(2*pi*40*t) + Au*sin(2*pi*60*t);
[ts,xsp9,ys] = sim('modelP9',TIME_SIM,[],[t' ue']); 
xsrp9 = xsp9(1/STEP_SIZE:length(t),:); %xsr = xs; % xsr = xs(1/STEP_SIZE:length(t),:);

figure
grid on, axP = axes; set(axP, 'FontSize', 14)
subplot(411), plot(t(1:20000),ue(1:20000))
title('input $u_e$ and states responses in time domain', 'FontSize', 14, 'Interpreter','Latex')
xlabel('Time [s]', 'FontSize', 14, 'Interpreter','Latex')
ylabel('$u_e$ [V]', 'FontSize', 14, 'Interpreter','Latex')
subplot(412), plot(t(1:20000),xsp9(1:20000,1)); hold on; plot(t(1:20000),xsp3(1:20000,1),'--');
xlabel('Time [s]', 'FontSize', 14, 'Interpreter','Latex')
ylabel('$x$ [m]', 'FontSize', 14, 'Interpreter','Latex')
l = legend('linear model', 'nonlinear model');
set(l, 'Interpreter','Latex');
subplot(413), plot(t(1:20000),xsp9(1:20000,2)); hold on; plot(t(1:20000),xsp3(1:20000,2),'--');
xlabel('Time [s]', 'FontSize', 14, 'Interpreter','Latex')
ylabel('$\dot{x}$ [m/s]', 'FontSize', 14, 'Interpreter','Latex')
l = legend('linear model', 'nonlinear model');
set(l, 'Interpreter','Latex');
subplot(414), plot(t(1:20000),xsp9(1:20000,3)); hold on; plot(t(1:20000),xsp3(1:20000,3),'--');
xlabel('Time [s]', 'FontSize', 14, 'Interpreter','Latex')
ylabel('$i$ [A]', 'FontSize', 14, 'Interpreter','Latex')
l = legend('linear model', 'nonlinear model');
set(l, 'Interpreter','Latex');

% freq domain

Pxxp9=[];
fp9=[];
for i = 1:1:3
    [Pxx_tmp, f_tmp]=power_spectral_density(xsrp9(:,i), Fs);
    Pxxp9=[Pxxp9, Pxx_tmp];
    fp9=[fp9, f_tmp'];
end
[UE, fe]=power_spectral_density(uer, Fs);

Pxxp9_db=10*log10(Pxxp9);
UE_db=10*log10(UE);

figure
grid on, axP = axes; set(axP, 'FontSize', 14)
subplot(411), plot(fe,UE_db)
title('Periodogram Using FFT')
axis([fc-10 5*fc -inf inf]);
xlabel('frequency [Hz]', 'FontSize', 14, 'Interpreter','Latex')
ylabel('$P_{u_e}$ [dB/Hz]', 'FontSize', 14, 'Interpreter','Latex')
subplot(412), plot(fp9(:,1),Pxxp9_db(:,1));% hold on; plot(fp9(:,1),Pxxp3_db(:,1))
axis([fc-10 5*fc -inf inf]);
xlabel('frequency [Hz]', 'FontSize', 14, 'Interpreter','Latex')
ylabel('$P_{x}$ [dB/Hz]', 'FontSize', 14, 'Interpreter','Latex')
%l = legend('linear model', 'nonlinear model');
%set(l, 'Interpreter','Latex');
subplot(413), plot(fp9(:,2),Pxxp9_db(:,2));% hold on; plot(fp9(:,2),Pxxp3_db(:,2))
axis([fc-10 5*fc -inf inf]);
xlabel('frequency [Hz]', 'FontSize', 14, 'Interpreter','Latex')
ylabel('$P_{\dot{x}}$ [dB/Hz]', 'FontSize', 14, 'Interpreter','Latex')
%l = legend('linear model', 'nonlinear model');
%set(l, 'Interpreter','Latex');
subplot(414), plot(fp9(:,3),Pxxp9_db(:,3));% hold on; plot(fp9(:,3),Pxxp3_db(:,3))
axis([fc-10 5*fc -inf inf]);
xlabel('frequency [Hz]', 'FontSize', 14, 'Interpreter','Latex')
ylabel('$P_i$ [dB/Hz]', 'FontSize', 14, 'Interpreter','Latex')
%l = legend('linear model', 'nonlinear model');
%set(l, 'Interpreter','Latex');

%% PB10
C=[0 0 1];
sys=ss(A,B,C,[0]);
w = [2:2:6]*pi*fc;
[MAG,PHASE] = bode(sys,w);
MAGp10=[MAG(1) MAG(2) MAG(3)];
NoiseAmp = Amplitudep4(3,1:3)./MAGp10;

Noise1 = NoiseAmp(2)*sin(4*pi*fc*t);
Noise2 = NoiseAmp(3)*sin(6*pi*fc*t);
Noise = [Noise1; Noise2];

Bd=[B B];
[ts,xsp10,ys] = sim('modelP11',TIME_SIM,[],[t' ue'], [t' Noise1'], [t' Noise2']); 
xsrp10 = xsp10(1/STEP_SIZE:length(t),:); 

%time domain
grid on, axP = axes; set(axP, 'FontSize', 14)
subplot(411), plot(t(1:20000),ue(1:20000))
title('input $u_e$ and states responses in time domain', 'FontSize', 14, 'Interpreter','Latex')
xlabel('Time [s]', 'FontSize', 14, 'Interpreter','Latex')
ylabel('$u_e$ [V]', 'FontSize', 14, 'Interpreter','Latex')
subplot(412), plot(t(1:20000),xsp10(1:20000,1)); hold on; plot(t(1:20000),xsp3(1:20000,1),'--');
xlabel('Time [s]', 'FontSize', 14, 'Interpreter','Latex')
ylabel('$x$ [m]', 'FontSize', 14, 'Interpreter','Latex')
l = legend('linear model', 'nonlinear model');
set(l, 'Interpreter','Latex');
subplot(413), plot(t(1:20000),xsp10(1:20000,2)); hold on; plot(t(1:20000),xsp3(1:20000,2),'--');
xlabel('Time [s]', 'FontSize', 14, 'Interpreter','Latex')
ylabel('$\dot{x}$ [m/s]', 'FontSize', 14, 'Interpreter','Latex')
l = legend('linear model', 'nonlinear model');
set(l, 'Interpreter','Latex');
subplot(414); 
plot(t(1:20000),xsp10(1:20000,3)); hold on; plot(t(1:20000),xsp3(1:20000,3),'--');
xlabel('Time [s]', 'FontSize', 14, 'Interpreter','Latex')
ylabel('$i$ [A]', 'FontSize', 14, 'Interpreter','Latex')
l = legend('linear model', 'nonlinear model');
set(l, 'Interpreter','Latex');


% freq domain
Pxxp10=[];
fp10=[];
for i = 1:1:3
[Pxx_tmp, f_tmp]=power_spectral_density(xsrp10(:,i), Fs);
Pxxp10=[Pxxp10, Pxx_tmp];
fp10=[fp10, f_tmp'];
end
[UE, fe]=power_spectral_density(uer, Fs);
Pxxp10_db=10*log10(Pxxp10);
UE_db=10*log10(UE);
figure
grid on, axP = axes; set(axP, 'FontSize', 14)
subplot(411), plot(fe,UE_db)
title('Periodogram Using FFT')
axis([fc-10 5*fc -inf inf]);
xlabel('frequency [Hz]', 'FontSize', 14, 'Interpreter','Latex')
ylabel('$P_{u_e}$ [dB/Hz]', 'FontSize', 14, 'Interpreter','Latex')
subplot(412), plot(fp10(:,1),Pxxp10_db(:,1)); hold on; plot(fp3(:,1),Pxxp3_db(:,1));
l = legend('linear model', 'nonlinear model');
set(l, 'Interpreter','Latex');
axis([fc-10 5*fc -inf inf]);
xlabel('frequency [Hz]', 'FontSize', 14, 'Interpreter','Latex')
ylabel('$P_{x}$ [dB/Hz]', 'FontSize', 14, 'Interpreter','Latex')
subplot(413), plot(fp10(:,2),Pxxp10_db(:,2)); hold on; plot(fp3(:,2),Pxxp3_db(:,2));
l = legend('linear model', 'nonlinear model');
set(l, 'Interpreter','Latex');
axis([fc-10 5*fc -inf inf]);
xlabel('frequency [Hz]', 'FontSize', 14, 'Interpreter','Latex')
ylabel('$P_{\dot{x}}$ [dB/Hz]', 'FontSize', 14, 'Interpreter','Latex')
subplot(414); plot(fp10(:,3),Pxxp10_db(:,3)); hold on; plot(fp3(:,3),Pxxp3_db(:,3),'--');
l = legend('linear model', 'nonlinear model');
set(l, 'Interpreter','Latex');
axis([fc-10 5*fc -inf inf]);
xlabel('frequency [Hz]', 'FontSize', 14, 'Interpreter','Latex')
ylabel('$P_{i}$ [dB/Hz]', 'FontSize', 14, 'Interpreter','Latex')
Pxxp10_db(40*TIME_SIM+1, 3), Pxxp10_db(60*TIME_SIM+1, 3);

amplitude(xsp10(:,3), TIME_SIM, 40);
amplitude(xsp10(:,3), TIME_SIM, 60);

%%
%Pb12
Ts=0.0002; %[s]
[F,G]=c2d(A, B, Ts);
[F, Gd] = c2d(A, Bd, Ts);
lambdad=eig(F);

%% Pb13
k=0:TIME_SIM/STEP_SIZE;
ued=Au*sin(2*pi*fc*k*Ts);
Noise1d = NoiseAmp(2)*sin(4*pi*fc*k*Ts);
Noise2d = NoiseAmp(3)*sin(6*pi*fc*k*Ts);

[ts,xsp13,ys] = sim('modelP13',TIME_SIM,[],[(k*Ts)' ued'], [(k*Ts)' Noise1d'], [(k*Ts)' Noise2d']);
uedr = uoutP13.Data(1/Ts:end);   
xsrp13 = xoutP13.Data(1/Ts:end,:);
kr = k(1/Ts:length(k));
%time domain
figure
grid on, axP = axes; set(axP, 'FontSize', 14)
subplot(411), stairs(k(1:2/Ts)*Ts,uoutP13.Data(1:2/Ts)); hold on; plot(t(1:2/STEP_SIZE),ue(1:2/STEP_SIZE),'--');
title('input $u_e$ and states responses in time domain', 'FontSize', 14, 'Interpreter','Latex')
xlabel('Time [s]', 'FontSize', 14, 'Interpreter','Latex')
ylabel('$u_e$ [V]', 'FontSize', 14, 'Interpreter','Latex')
l = legend('Discrete Time model', 'Continuous Time Model');
set(l, 'Interpreter','Latex');
subplot(412), stairs(k(1:2/Ts)*Ts,xoutP13.Data((1:2/Ts),1)); hold on; plot(t(1:2/STEP_SIZE),xsp10((1:2/STEP_SIZE),1),'--');
xlabel('Time [s]', 'FontSize', 14, 'Interpreter','Latex')
ylabel('$x$ [m]', 'FontSize', 14, 'Interpreter','Latex')
l = legend('Discrete Time model', 'Continuous Time Model');
set(l, 'Interpreter','Latex');
subplot(413), stairs(k(1:2/Ts)*Ts,xoutP13.Data((1:2/Ts),2)); hold on; plot(t(1:2/STEP_SIZE),xsp10((1:2/STEP_SIZE),2),'--');
xlabel('Time [s]', 'FontSize', 14, 'Interpreter','Latex')
ylabel('$\dot{x}$ [m/s]', 'FontSize', 14, 'Interpreter','Latex')
l = legend('Discrete Time model', 'Continuous Time Model');
set(l, 'Interpreter','Latex');
subplot(414), stairs(k(1:2/Ts)*Ts,xoutP13.Data((1:2/Ts),3)); hold on; plot(t(1:2/STEP_SIZE),xsp10((1:2/STEP_SIZE),3),'--');
xlabel('Time [s]', 'FontSize', 14, 'Interpreter','Latex')
ylabel('$i$ [A]', 'FontSize', 14, 'Interpreter','Latex')
l = legend('Discrete Time model', 'Continuous Time Model');
set(l, 'Interpreter','Latex');

%freq domain
Pxxp13=[];
fp13=[];
Fsd = 1/Ts;
for i = 1:1:3
[Pxx_tmp, f_tmp]=power_spectral_density(xsrp13(:,i), Fsd);
Pxxp13=[Pxxp13, Pxx_tmp];
fp13=[fp13, f_tmp'];
end
[UEd, fed]=power_spectral_density(uedr, Fsd);
Pxxp13_db=10*log10(Pxxp13);
UE_dbd=10*log10(UEd);

%freq domain
figure
grid on, axP = axes; set(axP, 'FontSize', 14)
subplot(411), plot(fed,UE_dbd); hold on; plot(fe,UE_db,'--');
title('Periodogram Using FFT')
axis([fc-10 5*fc -inf inf]);
xlabel('frequency [Hz]', 'FontSize', 14, 'Interpreter','Latex')
ylabel('$P_{u_e}$ [dB/Hz]', 'FontSize', 14, 'Interpreter','Latex')
l = legend('Discrete Time model', 'Continuous Time Model');
set(l, 'Interpreter','Latex');
subplot(412), plot(fp13(:,1),Pxxp13_db(:,1)); hold on; plot(fp10(:,1),Pxxp10_db(:,1),'--');
axis([fc-10 5*fc -inf inf]);
xlabel('frequency [Hz]', 'FontSize', 14, 'Interpreter','Latex')
ylabel('$P_{x}$ [dB/Hz]', 'FontSize', 14, 'Interpreter','Latex')
l = legend('Discrete Time model', 'Continuous Time Model');
set(l, 'Interpreter','Latex');
subplot(413), plot(fp13(:,2),Pxxp13_db(:,2)); hold on; plot(fp10(:,2),Pxxp10_db(:,2),'--');
axis([fc-10 5*fc -inf inf]);
xlabel('frequency [Hz]', 'FontSize', 14, 'Interpreter','Latex')
ylabel('$P_{\dot{x}}$ [dB/Hz]', 'FontSize', 14, 'Interpreter','Latex')
l = legend('Discrete Time model', 'Continuous Time Model');
set(l, 'Interpreter','Latex');
subplot(414), plot(fp13(:,3),Pxxp13_db(:,3)); hold on; plot(fp10(:,3),Pxxp10_db(:,3),'--');
axis([fc-10 5*fc -inf inf]);
xlabel('frequency [Hz]', 'FontSize', 14, 'Interpreter','Latex')
ylabel('$P_{i}$ [dB/Hz]', 'FontSize', 14, 'Interpreter','Latex')
l = legend('Discrete Time model', 'Continuous Time Model');
set(l, 'Interpreter','Latex');


%% Pb14
Mc=[G F*G F^2*G];
rank(Mc); % = 3, controllable

%% Pb15
C=[1 0 0];
%continous to discrete
Pconti= [-4000, -1200-1000i, -1200+1000i];%[-2000, -1200-1000i, -1200+1000i]; %[-900, -700-400i, -700+400i];
% Kcon=acker(A,B,Pconti);
% Ak=A-B*Kcon;
% sys=ss(A,B,C,[]);
% sys_k=ss(Ak,B,C,[]);
% w = 2*pi*fc;
% figure;
% bode(sys_k, [0:0.1:10000]);
% [MAG2,PHASE] = bode(sys,w);
% Pdisc = exp(Pconti*Ts);
% 
% K=acker(F,G,Pdisc);
% Fk=F-G*K;
% sys_d=ss(F,G,C,[]);
Pconti = Pconti*1.4;
% discrete
[z p kss] = ss2zp(F,G,C,0);
%Pdisc = [p(1)*0.9 p(2)*0.9 p(3)*0.9];
Pdisc = exp(Pconti*Ts)
%Pdisc = [0.6703 + 0.0000i, 0.8825 - 0.0885i, 0.8825 + 0.0885i];
K=acker(F,G,Pdisc);
Fk=F-G*K;
sys_kd=ss(Fk,G,C,0,Ts);
bode(sys_kd);
tau1 = 1/Pconti(1)
ome1 = sqrt(imag(Pconti(2))^2 + real(Pconti(2))^2)
dzeta1 = -real(Pconti(2))/ome1
ome1/(2*pi)

%%
Ax=0.0008; %[m]
%N=inv(C*inv(-Ak)*B);
N = inv(C*inv(eye(3)-Fk)*G);
%Av=Ax*N; %8.3749; %Ax/MAG2;
xref=Ax*sin(2*pi*fc*k*Ts);
[ts,xsp15,ys] = sim('modelP15',TIME_SIM,[],[(k*Ts)' xref'], [(k*Ts)' Noise1d'], [(k*Ts)' Noise2d']);
xrefr = xref(1/Ts:end);   
xsrp15 = xsp15(1/Ts:end,:);

%time domain
figure
grid on, axP = axes; set(axP, 'FontSize', 14)
subplot(411),plot((k(2/Ts:2.2/Ts)*Ts),xoutP15.Data((2/Ts:2.2/Ts),1)); hold on; plot((k(2/Ts:2.2/Ts)*Ts),xref(2/Ts:2.2/Ts),'--')
l = legend('state $x$', '$x_{ref}$');
set(l, 'Interpreter','Latex');
xlabel('Time [s]', 'FontSize', 14, 'Interpreter','Latex')
ylabel('$x$ [m]', 'FontSize', 14, 'Interpreter','Latex')

figure;
grid on, axP = axes; set(axP, 'FontSize', 14)
subplot(411), plot(k(1:2/Ts)*Ts,ucontP15.Data(1:2/Ts))
title('control input $u_e$', 'FontSize', 14, 'Interpreter','Latex')
xlabel('Time [s]', 'FontSize', 14, 'Interpreter','Latex')
ylabel('$u_e$', 'FontSize', 14, 'Interpreter','Latex')


%% Pb16
Ax=0.001;
xref_c1 = Ax*sin(2*pi*fc*t);
xref_c0 = 0.0005*sin(2*pi*20*t);
Noise10 = 0*t;
Noise20 = 0*t;
[ts,xsp1601,ys] = sim('modelP16',TIME_SIM,[],[t' xref_c1'], [t' Noise10'], [t' Noise20']);
xref_cr = xref_c1(1/STEP_SIZE:length(t));   
xsrp1601 = xsp1601(1/STEP_SIZE:length(t),:);

tmp2=ucontP16.Data;

%xsp160=xsp1601;
% Ax=0.00075;
% xref_c2 = Ax*sin(2*pi*fc*t);
% [ts,xsp1602,ys] = sim('modelP16',TIME_SIM,[],[t' xref_c2'], [t' Noise10'], [t' Noise20']);
% Ax=0.001;
% xref_c3 = Ax*sin(2*pi*fc*t);
% [ts,xsp1603,ys] = sim('modelP16',TIME_SIM,[],[t' xref_c3'], [t' Noise10'], [t' Noise20']);
%time domain
figure
grid on, axP = axes; set(axP, 'FontSize', 14)
% subplot(411), plot(t(1:20000),xref_c(1:20000))
% title('input $u_e$ and states responses in time domain', 'FontSize', 14, 'Interpreter','Latex')
% xlabel('Time [s]', 'FontSize', 14, 'Interpreter','Latex')
% ylabel('$u_e$', 'FontSize', 14, 'Interpreter','Latex')
% subplot(412), plot(t(1:20000),xsp160(1:20000,1)); hold on; plot(t(1:20000),xref_c(1:20000));
% l = legend('state $x$', '$x_{ref}$');
% set(l, 'Interpreter','Latex');
% xlabel('Time [s]', 'FontSize', 14, 'Interpreter','Latex')
% ylabel('$x$ [m]', 'FontSize', 14, 'Interpreter','Latex')
% subplot(413), plot(t(1:20000),xsp160(1:20000,2))
% xlabel('Time [s]', 'FontSize', 14, 'Interpreter','Latex')
% ylabel('$\dot{x}$ [m/s]', 'FontSize', 14, 'Interpreter','Latex')
subplot(211), plot(t(2/STEP_SIZE:2.1/STEP_SIZE),xsp160(2/STEP_SIZE:2.1/STEP_SIZE,1)); hold on; plot(t(2/STEP_SIZE:2.1/STEP_SIZE),xref_c0(2/STEP_SIZE:2.1/STEP_SIZE),'--');
l = legend('state $x$', '$x_{ref}$');
set(l, 'Interpreter','Latex');
xlabel('Time [s]', 'FontSize', 14, 'Interpreter','Latex')
ylabel('$x$ [m]', 'FontSize', 14, 'Interpreter','Latex')
subplot(212), plot(t(2/STEP_SIZE:2.1/STEP_SIZE),xsp1601(2/STEP_SIZE:2.1/STEP_SIZE,1)); hold on; plot(t(2/STEP_SIZE:2.1/STEP_SIZE),xref_c1(2/STEP_SIZE:2.1/STEP_SIZE),'--');
l = legend('state $x$', '$x_{ref}$');
set(l, 'Interpreter','Latex');
xlabel('Time [s]', 'FontSize', 14, 'Interpreter','Latex')
ylabel('$x$ [m]', 'FontSize', 14, 'Interpreter','Latex')
% subplot(313), plot(t(1:20000),xsp1603(1:20000,1)); hold on; plot(t(1:20000),xref_c3(1:20000),'--');
% l = legend('state $x$', '$x_{ref}$');
% set(l, 'Interpreter','Latex');
% xlabel('Time [s]', 'FontSize', 14, 'Interpreter','Latex')
% ylabel('$x$ [m]', 'FontSize', 14, 'Interpreter','Latex')


%non zero noise
Ax=0.001;
%xsp91=xsp9;
xref_c1 = Ax*sin(2*pi*fc*t);
[ts,xsp161,ys] = sim('modelP16',TIME_SIM,[],[t' xsp9(:,1)], [t' Noise1'], [t' Noise2']);
xsrp161 = xsp161(1/STEP_SIZE:length(t),:);
%xsp161tmp=xsp161;
% Ax=0.00075;
% xref_c2 = Ax*sin(2*pi*fc*t);
% [ts,xsp162,ys] = sim('modelP16',TIME_SIM,[],[t' xref_c2'], [t' Noise1'], [t' Noise2']);
% Ax=0.001;
% xref_c3 = Ax*sin(2*pi*fc*t);
% [ts,xsp163,ys] = sim('modelP16',TIME_SIM,[],[t' xref_c3'], [t' Noise1'], [t' Noise2']);
%time domain
figure
grid on, axP = axes; set(axP, 'FontSize', 14)
% subplot(411), plot(t(1:20000),xref_c(1:20000))
% title('input $u_e$ and states responses in time domain', 'FontSize', 14, 'Interpreter','Latex')
% xlabel('Time [s]', 'FontSize', 14, 'Interpreter','Latex')
% ylabel('$u_e$', 'FontSize', 14, 'Interpreter','Latex')
subplot(211), plot(t(2/STEP_SIZE:2.1/STEP_SIZE),xsp161tmp(2/STEP_SIZE:2.1/STEP_SIZE,1)); hold on; plot(t(2/STEP_SIZE:2.1/STEP_SIZE),xsp91(2/STEP_SIZE:2.1/STEP_SIZE,1)','--');
l = legend('state $x$', '$x_{ref}$');
set(l, 'Interpreter','Latex');
xlabel('Time [s]', 'FontSize', 14, 'Interpreter','Latex')
ylabel('$x$ [m]', 'FontSize', 14, 'Interpreter','Latex')
subplot(212), plot(t(2/STEP_SIZE:2.1/STEP_SIZE),xsp161(2/STEP_SIZE:2.1/STEP_SIZE,1)); hold on; plot(t(2/STEP_SIZE:2.1/STEP_SIZE),xsp9(2/STEP_SIZE:2.1/STEP_SIZE,1)','--');
l = legend('state $x$', '$x_{ref}$');
set(l, 'Interpreter','Latex');
xlabel('Time [s]', 'FontSize', 14, 'Interpreter','Latex')
ylabel('$x$ [m]', 'FontSize', 14, 'Interpreter','Latex')
% subplot(312), plot(t(1:20000),xsp162(1:20000,1)); hold on; plot(t(1:20000),xref_c2(1:20000),'--');
% l = legend('state $x$', '$x_{ref}$');
% set(l, 'Interpreter','Latex');
% xlabel('Time [s]', 'FontSize', 14, 'Interpreter','Latex')
% ylabel('$x$ [m]', 'FontSize', 14, 'Interpreter','Latex')
% subplot(313), plot(t(1:20000),xsp163(1:20000,1)); hold on; plot(t(1:20000),xref_c3(1:20000),'--');
% l = legend('state $x$', '$x_{ref}$');
% set(l, 'Interpreter','Latex');
% xlabel('Time [s]', 'FontSize', 14, 'Interpreter','Latex')
% ylabel('$x$ [m]', 'FontSize', 14, 'Interpreter','Latex')

% subplot(413), plot(t(1:20000),xsp16(1:20000,2))
% xlabel('Time [s]', 'FontSize', 14, 'Interpreter','Latex')
% ylabel('$\dot{x}$ [m/s]', 'FontSize', 14, 'Interpreter','Latex')
% subplot(414), plot(t(1:20000),xsp16(1:20000,3))
% xlabel('Time [s]', 'FontSize', 14, 'Interpreter','Latex')
% ylabel('$i$ [A]', 'FontSize', 14, 'Interpreter','Latex')

% Pxxp1601=[];
% fp1601=[];
% for i = 1:1:3
% [Pxx_tmp, f_tmp]=power_spectral_density(xsrp1601(:,i), Fs);
% Pxxp1601=[Pxxp1601, Pxx_tmp];
% fp1601=[fp1601, f_tmp'];
% end
% [XREF, fe]=power_spectral_density(xref_cr, Fs);
% Pxxp1601_db=10*log10(Pxxp1601);
% XREF_db=10*log10(XREF);
% 
% Pxxp161=[];
% fp161=[];
% for i = 1:1:3
% [Pxx_tmp, f_tmp]=power_spectral_density(xsrp161(:,i), Fs);
% Pxxp161=[Pxxp161, Pxx_tmp];
% fp161=[fp161, f_tmp'];
% end
% [XREF, fe]=power_spectral_density(xref_cr, Fs);
% Pxxp161_db=10*log10(Pxxp161);
% XREF_db=10*log10(XREF);
% 
% figure
% grid on, axP = axes; set(axP, 'FontSize', 14)
% subplot(211), plot(fp1601(:,3),Pxxp1601_db(:,3)); hold on; plot(fp3(:,3),Pxxp3_db(:,3));
% title('Periodogram Using FFT')
% l = legend('linear model', 'nonlinear model');
% set(l, 'Interpreter','Latex');
% axis([fc-10 5*fc -inf inf]);
% xlabel('frequency [Hz]', 'FontSize', 14, 'Interpreter','Latex')
% ylabel('PSD of $i$ [dB/Hz]', 'FontSize', 14, 'Interpreter','Latex')
% subplot(212), plot(fp161(:,3),Pxxp161_db(:,3)); hold on; plot(fp3(:,3),Pxxp3_db(:,3));
% title('Periodogram Using FFT')
% l = legend('linear model', 'nonlinear model');
% set(l, 'Interpreter','Latex');
% axis([fc-10 5*fc -inf inf]);
% xlabel('frequency [Hz]', 'FontSize', 14, 'Interpreter','Latex')
% ylabel('PSD of $i$ [dB/Hz]', 'FontSize', 14, 'Interpreter','Latex')




figure;
grid on
subplot(211); plot(t(2/STEP_SIZE:2.1/STEP_SIZE),tmp2(2/STEP_SIZE:2.1/STEP_SIZE))
title('control input $u_e$', 'Interpreter','Latex')
xlabel('Time [s]', 'Interpreter','Latex')
ylabel('$u_e$', 'Interpreter','Latex')
subplot(212); plot(t(2/STEP_SIZE:2.1/STEP_SIZE),ucontP16.Data(2/STEP_SIZE:2.1/STEP_SIZE))
title('control input $u_e$', 'Interpreter','Latex')
xlabel('Time [s]', 'Interpreter','Latex')
ylabel('$u_e$', 'Interpreter','Latex')

%% P17
Ax=0.001;
xref_c = Ax*sin(2*pi*fc*t);

[ts,xsp17,ys] = sim('modelP17',TIME_SIM,[],[t' xref_c']);
xsrp17 = xsp17(1/STEP_SIZE:length(t),:);
xref_cr = xref_c(1/STEP_SIZE:length(t));  

%time domain
figure
grid on, axP = axes; set(axP, 'FontSize', 14)
subplot(411), plot(t(1:20000),xsp17(1:20000,1)); hold on; plot(t(1:20000),xref_c(1:20000),'--');
l = legend('state $x$', '$x_{ref}$');
set(l, 'Interpreter','Latex');
xlabel('Time [s]', 'FontSize', 14, 'Interpreter','Latex')
ylabel('$x$ [m]', 'FontSize', 14, 'Interpreter','Latex')
subplot(412), plot(t(1:20000),xsp17(1:20000,2))
xlabel('Time [s]', 'FontSize', 14, 'Interpreter','Latex')
ylabel('$\dot{x}$ [m/s]', 'FontSize', 14, 'Interpreter','Latex')
subplot(413), plot(t(1:20000),xsp17(1:20000,3))
xlabel('Time [s]', 'FontSize', 14, 'Interpreter','Latex')
ylabel('$i$ [A]', 'FontSize', 14, 'Interpreter','Latex')

subplot(414), plot(t(1:20000),ucontP17.Data(1:20000))
title('control input $u_e$', 'FontSize', 14, 'Interpreter','Latex')
xlabel('Time [s]', 'FontSize', 14, 'Interpreter','Latex')
ylabel('$u_e$ [V]', 'FontSize', 14, 'Interpreter','Latex')

%freq domain
Pxxp17=[];
fp17=[];
for i = 1:1:3
[Pxx_tmp, f_tmp]=power_spectral_density(xsrp17(:,i), Fs);
Pxxp17=[Pxxp17, Pxx_tmp];
fp17=[fp17, f_tmp'];
end
[XREF, fe]=power_spectral_density(xref_cr, Fs);
Pxxp17_db=10*log10(Pxxp17);
XREF_db=10*log10(XREF);

figure
grid on, axP = axes; set(axP, 'FontSize', 14)
subplot(411), plot(fe,XREF_db)
title('Periodogram Using FFT')
axis([fc-10 5*fc -inf inf]);
xlabel('frequency [Hz]', 'FontSize', 14, 'Interpreter','Latex')
ylabel('$P_{u_e}$ [dB/Hz]', 'FontSize', 14, 'Interpreter','Latex')
subplot(412), plot(fp17(:,1),Pxxp17_db(:,1)); hold on; plot(fp3(:,1),Pxxp3_db(:,1));
l = legend('nonlinear model with controller', 'nonlinear model');
set(l, 'Interpreter','Latex');
axis([fc-10 5*fc -inf inf]);
xlabel('frequency [Hz]', 'FontSize', 14, 'Interpreter','Latex')
ylabel('$P_x$ [dB/Hz]', 'FontSize', 14, 'Interpreter','Latex')
subplot(413), plot(fp17(:,2),Pxxp17_db(:,2)); hold on; plot(fp3(:,2),Pxxp3_db(:,2));
l = legend('nonlinear model with controller', 'nonlinear model');
set(l, 'Interpreter','Latex');
axis([fc-10 5*fc -inf inf]);
xlabel('frequency [Hz]', 'FontSize', 14, 'Interpreter','Latex')
ylabel('$P_{\dot{x}}$ [dB/Hz]', 'FontSize', 14, 'Interpreter','Latex')
subplot(414), plot(fp17(:,3),Pxxp17_db(:,3)); hold on; plot(fp3(:,3),Pxxp3_db(:,3));
l = legend('nonlinear model with controller', 'nonlinear model');
set(l, 'Interpreter','Latex');
axis([fc-10 5*fc -inf inf]);
xlabel('frequency [Hz]', 'FontSize', 14, 'Interpreter','Latex')
ylabel('$P_i$ [dB/Hz]', 'FontSize', 14, 'Interpreter','Latex')

Amplitudep17=[];
for j=1:3
    for i=1:6
       Amplitudep17(j,i)=amplituder(xsrp17(:,j), TIME_SIM, fc*i);
    end
end

numberHarmo = 6;
THDp17 = THD(Amplitudep17(2,:),numberHarmo);
d2_control = Amplitudep17(2,2)/sqrt(Amplitudep17(2,1)^2+Amplitudep17(2,2)^2)*100;
d3_control = Amplitudep17(2,3)/sqrt(Amplitudep17(2,1)^2+Amplitudep17(2,3)^2)*100;

%% Pb 18
% AxPb18 = 0.0005:0.00025:0.001; % m
% fcPb18 = 20:5:200; % Hz
% d2Pb18 = zeros(length(AxPb18),length(fcPb18));
% d3Pb18 = zeros(length(AxPb18),length(fcPb18));
% 
% for k=1:length(AxPb18)
%     for l=1:length(fcPb18)
%         xrefPb18 = AxPb18(k)*sin(2*pi*fcPb18(l)*t);
%         [tsPb18,xsPb18,ysPb18] = sim('modelP17',TIME_SIM,[],[t' xrefPb18']);
%         xsrPb18 = xsPb18(1/STEP_SIZE:length(t),:); 
%         
%         AmplitudePb18=[];
%         for j=1:3
%             for i=1:3
%                AmplitudePb18(j,i)=amplituder(xsrPb18(:,j), TIME_SIM, fcPb18(l)*i);
%             end
%         end
%         d2Pb18(k,l) = AmplitudePb18(2,2)/sqrt(AmplitudePb18(2,1)^2+AmplitudePb18(2,2)^2)*100;
%         d3Pb18(k,l) = AmplitudePb18(2,3)/sqrt(AmplitudePb18(2,1)^2+AmplitudePb18(2,3)^2)*100;
%     end
% end
% 
% figure
% grid on, axP = axes; set(axP, 'FontSize', 14)
% subplot(311)
% plot(fcPb18, d2Pb18(1,:))
% hold on
% plot(fcPb18, d3Pb18(1,:))
% xlabel('Frequency [Hz]', 'FontSize', 14, 'Interpreter','Latex')
% ylabel('$A_x = 0.0005$ m', 'FontSize', 14, 'Interpreter','Latex')
% l=legend('$d_2 [\%]$','$d_3 [\%]$');
% set(l, 'FontSize',14, 'Interpreter','Latex')
% subplot(312)
% plot(fcPb18, d2Pb18(2,:))
% hold on
% plot(fcPb18, d3Pb18(2,:))
% xlabel('Frequency [Hz]', 'FontSize', 14, 'Interpreter','Latex')
% ylabel('$A_x = 0.00075$ m', 'FontSize', 14, 'Interpreter','Latex')
% l=legend('$d_2 [\%]$','$d_3 [\%]$');
% set(l, 'FontSize',14, 'Interpreter','Latex')
% subplot(313)
% plot(fcPb18, d2Pb18(3,:))
% hold on
% plot(fcPb18, d3Pb18(3,:))
% xlabel('Frequency [Hz]', 'FontSize', 14, 'Interpreter','Latex')
% ylabel('$A_x = 0.001$ m', 'FontSize', 14, 'Interpreter','Latex')
% l=legend('$d_2 [\%]$','$d_3 [\%]$');
% set(l, 'FontSize',14, 'interpreter','Latex')


%% Pb 19
C_p19=[0 0 1];

%TODO
% Axw = zeros(3,4);
% Aw = zeros(4,4);
% 
% 
% %should be ok
% Ak_p19 = [A Axw;
%     zeros(4,3) Aw];
% Bk_p19 = [B;0;0;0;0];
% Ck_p19 = [C_p19 0 0 0 0];
% 
% 
% Mo_p19 = [Ck_p19
%     Ck_p19*Ak_p19
%     Ck_p19*Ak_p19^2
%     Ck_p19*Ak_p19^3
%     Ck_p19*Ak_p19^4
%     Ck_p19*Ak_p19^5
%     Ck_p19*Ak_p19^6];
% rank(Mo_p19); % != 7 but observable

syms s var;
w40 = 2*pi*2*fc;
[num40, den40] = numden(laplace(sin(w40*var),s));
[Aw40 Bw40 Cw40 Dw40] = tf2ss(sym2poly(num40),sym2poly(den40));
w60 = 2*pi*3*fc;
[num60, den60] = numden(laplace(sin(w60*var),s));
[Aw60 Bw60 Cw60 Dw60] = tf2ss(sym2poly(num60),sym2poly(den60));
Aw = [Aw40 zeros(2,2)
    zeros(2,2) Aw60];
Axw = [Bd(:,1)*Cw40 Bd(:,2)*Cw60];
AkP19 = [A Axw;
    zeros(4,3) Aw];
lambdaAw = eig(Aw);
BkP19 = [B;zeros(4,1)];
CkP19 = [C_p19 0 0 0 0];
Akbis = [0 0;
        0 AkP19(2,2);
        AkP19(3,3) 0];
%sigmaVnd = [0.0002 0.0003 0.0004 0.0005];
sigmaVnd = [0 0 0 0];
Vnd = diag(sigmaVnd);
Bkn = [Akbis zeros(3,4); 
    zeros(4,2) eye(4)];

% observability
MoP19 = [CkP19
    CkP19*AkP19
    CkP19*AkP19^2
    CkP19*AkP19^3
    CkP19*AkP19^4
    CkP19*AkP19^5
    CkP19*AkP19^6];
rank(MoP19); % 5 != 7

%% Pb 20

sigmai = 0.1; % A
sigmav = 0.05; % m/s
sigman2 = 0.001; % A
beta = 125;
% sigma 2 = ajouter white noise a la sortie du syst = easy

% nx:
% cf solution GW15
% xa = [xP19
%        nx]
% with nxdot = -beta*nx+sqrt(2*beta)*[sigmai; sigmav]*whitenoise
% whitenoise ???

Aa = [AkP19 [Akbis; zeros(4,2)];
    zeros(2,7) -beta*eye(2)];

Ba = [BkP19 ; 0 ; 0];
Ca = [CkP19 0 0];
Bva = [zeros(3,6); [eye(4) zeros(4,2)]; [zeros(2,4) sqrt(2*beta)*diag([sigmai sigmav])]];

[Fa,Ga]=c2d(Aa, Ba, Ts);
[Fa, Gva] = c2d(Aa, Bva, Ts);
lambdap20=eig(Fa);

Mop20=[Ca
    Ca*Fa
    Ca*Fa^2
    Ca*Fa^3
    Ca*Fa^4
    Ca*Fa^5
    Ca*Fa^6
    Ca*Fa^7
    Ca*Fa^8];
rank(Mop20); % = 8 but observable

%% Pb 21
Q = Bva(4:end,:)*[Vnd zeros(4,2); zeros(2,4) diag([sigmai sigmav])]*Bva(4:end,:)'*Ts;
sigman2d = sigman2/Ts;
[M,P,Z,E] = dlqe(Fa,Gva,Ca,Q,sigman2d);

%% Pb 22
C = [0 0 1];
[ts,xsp22,ys] = sim('modelP21',TIME_SIM,[],[t' ue']); 

figure
grid on, axP = axes; set(axP, 'FontSize', 14)
subplot(411), plot(t(1:20000),xoutP22.Data(1:20000,1)); hold on; plot(k(1:20000)*Ts,xhatoutP22.Data(1:20000),'--');
l = legend('state $x$', '$x_{hat}$');
set(l, 'Interpreter','Latex');
xlabel('Time [s]', 'FontSize', 14, 'Interpreter','Latex')
ylabel('$x$ [m]', 'FontSize', 14, 'Interpreter','Latex')