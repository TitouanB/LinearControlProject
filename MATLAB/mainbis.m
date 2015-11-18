clear all; close all; clc;

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

[ts,xs,ys] = sim('nonLinearModel2',TIME_SIM,[],[t' ue']);

%%
figure
grid on, axP = axes; set(axP, 'FontSize', 14)
subplot(411), plot(t,ue)
title('input $u_e$ and states responses in time domain', 'FontSize', 14, 'Interpreter','Latex')
xlabel('Time [s]', 'FontSize', 14, 'Interpreter','Latex')
ylabel('$u_e$', 'FontSize', 14, 'Interpreter','Latex')
subplot(412), plot(t,xs(:,1))
xlabel('Time [s]', 'FontSize', 14, 'Interpreter','Latex')
ylabel('$x$ [m]', 'FontSize', 14, 'Interpreter','Latex')
subplot(413), plot(t,xs(:,2))
xlabel('Time [s]', 'FontSize', 14, 'Interpreter','Latex')
ylabel('$\dot{x}$ [m/s]', 'FontSize', 14, 'Interpreter','Latex')
subplot(414), plot(t,xs(:,3))
xlabel('Time [s]', 'FontSize', 14, 'Interpreter','Latex')
ylabel('$i$ [A]', 'FontSize', 14, 'Interpreter','Latex')

%%
Fs=1/STEP_SIZE;
Pxx=[];
f=[];
for i = 1:1:3
    [Pxx_tmp, f_tmp]=power_spectral_density(xs(:,i), Fs);
    Pxx=[Pxx, Pxx_tmp];
    f=[f, f_tmp'];
end
[UE, fe]=power_spectral_density(ue, Fs);

Pxx_db=10*log10(Pxx);
UE_db=10*log10(UE);

figure
grid on, axP = axes; set(axP, 'FontSize', 14)
subplot(411), plot(fe,UE_db)
title('Periodogram Using FFT')
axis([-1 100 -inf inf]);
xlabel('frequency [Hz]', 'FontSize', 14, 'Interpreter','Latex')
ylabel('PSD of $u_e$ [dB/Hz]', 'FontSize', 14, 'Interpreter','Latex')
subplot(412), plot(f(:,1),Pxx_db(:,1))
axis([-1 100 -inf inf]);
xlabel('frequency [Hz]', 'FontSize', 14, 'Interpreter','Latex')
ylabel('PSD of $x$ [dB/Hz]', 'FontSize', 14, 'Interpreter','Latex')
subplot(413), plot(f(:,2),Pxx_db(:,2))
axis([-1 100 -inf inf]);
xlabel('frequency [Hz]', 'FontSize', 14, 'Interpreter','Latex')
ylabel('PSD of $\dot{x}$ [dB/Hz]', 'FontSize', 14, 'Interpreter','Latex')
subplot(414), plot(f(:,3),Pxx_db(:,3))
axis([-1 100 -inf inf]);
xlabel('frequency [Hz]', 'FontSize', 14, 'Interpreter','Latex')
ylabel('PSD of $i$ [dB/Hz]', 'FontSize', 14, 'Interpreter','Latex')

%%
%PB4
%Aue=sqrt(Fs)*10.^(UE_db/20);

Aue = sqrt(10.^(UE_db(20*TIME_SIM+1)/10)*(20*length(ue)));
APxx = sqrt(10.^(Pxx_db/10)*Fs*length(xs));
%APxx=sqrt(Fs)*10.^(Pxx_db/20);



N = 6;
rue = thd(ue, Fs, N);
r1 = thd(xs(:,1), Fs, N);
r2 = thd(xs(:,2), Fs, N);
r3 = thd(xs(:,3), Fs, N);

% figure
% grid on, axP = axes; set(axP, 'FontSize', 14)
% subplot(411), plot(fe,Aue)
% axis([-1 100 -inf inf]);
% xlabel('frequency [Hz]', 'FontSize', 14, 'Interpreter','Latex')
% ylabel('PSD of $u_e$ [dB/Hz]', 'FontSize', 14, 'Interpreter','Latex')
% subplot(412), plot(f(:,1),APxx(:,1))
% axis([-1 100 -inf inf]);
% xlabel('frequency [Hz]', 'FontSize', 14, 'Interpreter','Latex')
% ylabel('PSD of $x$ [dB/Hz]', 'FontSize', 14, 'Interpreter','Latex')
% subplot(413), plot(f(:,2),APxx(:,2))
% axis([-1 100 -inf inf]);
% xlabel('frequency [Hz]', 'FontSize', 14, 'Interpreter','Latex')
% ylabel('PSD of $\dot{x}$ [dB/Hz]', 'FontSize', 14, 'Interpreter','Latex')
% subplot(414), plot(f(:,3),APxx(:,3))
% axis([-1 100 -inf inf]);
% xlabel('frequency [Hz]', 'FontSize', 14, 'Interpreter','Latex')
% ylabel('PSD of $i$ [dB/Hz]', 'FontSize', 14, 'Interpreter','Latex')

%% P6
[x,u,y,dx,options] = trim('nonLinearModel2',[0;0;0],[0],[],[],[1],[]);

[A,B,C,D] = linmod('nonLinearModelSubsystems',[0;0;0],0);