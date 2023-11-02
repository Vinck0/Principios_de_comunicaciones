%%
clear all
close all
clc
%% Parte 1
%Señal moduladora secuencia de Chirp
%----------------------------------------------------------------
%Señal up-chirp con frecuencia inicial f0 y frecuencia final.
%Se utiliza una resolución temporal para el barrido de 1 segundo.
fs_up=94.1*10^4*4;            %frecuencia de muestreo
t=0:1/fs_up:1;       %base temporal  para el ploteo de la señal
f0=1600;               %frecuencia inicial
f1=3000;       %frecuencia final

up=chirp_signal(t,f0,1,f1); %Función Chirp

figure()
subplot(2,1,1);
plot(t,up);
title('Up-Chirp Signal');
xlabel('Tiempo [s]');
ylabel('Amplitud[V]');

dt=1/fs_up; %Paso del Tiempo en segundos.
t=0:dt:1; %Base temporal de 1. para el ploteo de la señal.
N=length(t); %Largo del vector de Tiempo.

df=1/(N*dt); %Paso de la frecuencia en Hertz.
fshift=(-N/2:N/2-1)*df; %Vector de frecuencias centrado en 0.

subplot(2,1,2);
plot(fshift,fftshift(abs(fft(up)/fs_up))) %Gráfico de la señal modulante en el 
% dominio de la frecuencia, con la amplitud normalizada y centrada en 0
title('Up-Chirp Signal en el dominio de la frecuencia');
xlabel('frecuencia [Hz]');
ylabel('|F(Hz)|');
axis([1450 3050 0 inf])

%%
%----------------------------------------------------------------------------
%Señal down-chirp con frecuencia inicial f0 y f1 frecuencia final.
%Se utiliza una resolución temporal para el barrido de 1 segundo.
fs_down=94.1*10^4*4;            %frecuencia de muestreo
t=0:1/fs_down:1;       %base temporal  para el ploteo de la señal
f0=1500;               %frecuencia inicial
f1=300;         %frecuencia final

down=chirp_signal(t,f0,1,f1); %Función Chirp

figure()
subplot(2,1,1);
plot(t,down);
title('Down-Chirp Signal');
xlabel('Tiempo [s]');
ylabel('Amplitud[V]');

dt=1/fs_down; %Paso del Tiempo en segundos.
t=0:dt:1; %Base temporal de 0.2 segundos para el ploteo de la señal.
N=length(t); %Largo del vector de Tiempo.

df=1/(N*dt); %Paso de la frecuencia en Hertz.
fshift=(-N/2:N/2-1)*df; %Vector de frecuencias centrado en 0.
subplot(2,1,2);
plot(fshift,fftshift(abs(fft(down)/fs_down))) %Gráfico de la señal modulante en el 
% dominio de la frecuencia, con la amplitud normalizada y centrada en 0
title('Down-Chirp Signal en el dominio de la frecuencia');
xlabel('frecuencia [Hz]');
ylabel('|F(Hz)|');
axis([250 1550 0 inf])

%%
%---------------------------------------------------------------
%señal a modular = up + down chirp
senal = up + down;
figure()
subplot(2,1,1);
plot(t,senal);
title('Señal a ser Modulada');
xlabel('Tiempo [s]');
ylabel('Amplitud[V]');

dt=1/fs_down; %Paso del Tiempo en segundos.
t=0:dt:1; %Base temporal de 0.2 segundos para el ploteo de la señal.
N=length(t); %Largo del vector de Tiempo.

df=1/(N*dt); %Paso de la frecuencia en Hertz.
fshift=(-N/2:N/2-1)*df; %Vector de frecuencias centrado en 0.
subplot(2,1,2);
plot(fshift,fftshift(abs(fft(senal)/fs_up))) %Gráfico de la señal modulante en el 
% dominio de la frecuencia, con la amplitud normalizada y centrada en 0
title('Señal a ser modulada  en el dominio de la frecuencia');
xlabel('frecuencia [Hz]');
ylabel('|F(Hz)|');
axis([250 3100 0 inf])

%% Parte 2
fs=94.1*10^4*4;            %frecuencia de muestreo
t=0:1/fs:1;       %base temporal  para el ploteo de la señal
portadora = 1.5*cos(2*pi*94.1*10^4*t);

figure()
subplot(2,1,1);
plot(t,portadora) %Gráfico de la señal modulante en el 
% dominio de la frecuencia, con la amplitud normalizada y centrada en 0
title('Señal portadora');
xlabel('Tiempo [s]');
ylabel('Amplitud [V]');
axis([0 1 -2 2])

dt=1/fs; %Paso del Tiempo en segundos.
t=0:dt:1; %Base temporal de 0.2 segundos para el ploteo de la señal.
N=length(t); %Largo del vector de Tiempo.

df=1/(N*dt); %Paso de la frecuencia en Hertz.
fshift=(-N/2:N/2-1)*df; %Vector de frecuencias centrado en 0.
subplot(2,1,2);
plot(fshift,fftshift(abs(fft(portadora)/fs))) %Gráfico de la señal modulante en el 
% dominio de la frecuencia, con la amplitud normalizada y centrada en 0
title('Señal portadora en el dominio de la frecuencia');
xlabel('frecuencia [Hz]');
ylabel('|F(Hz)|');
axis([0 94.1*10^4*2 0 inf])


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%}


m1=100; %1000/94.1*10^6
m2=1000;
m3=10000;
m4=50000;
inte=dt*cumsum(senal);
%inte=trapz(0:1/fs_up:1,senal);
modulada1=1.5*cos(2*pi*94.1*10^4*t+2*pi*m1*inte);
modulada2=1.5*cos(2*pi*94.1*10^4*t+2*pi*m2*inte);
modulada3=1.5*cos(2*pi*94.1*10^4*t+2*pi*m3*inte);
modulada4=1.5*cos(2*pi*94.1*10^4*t+2*pi*m4*inte);

figure()
subplot(2,1,1);
plot(0:1/fs:1,modulada1); %Gráfico de la señal modulante en el 
% dominio de la frecuencia, con la amplitud normalizada y centrada en 0
title('Señal modulada f_d = 100[Hz]');
xlabel('Tiempo [s]');
ylabel('Amplitud [V]');
axis([0 inf -2 2])


subplot(2,1,2);
plot(fshift,fftshift(abs(fft(modulada1)/fs))) %Gráfico de la señal modulante en el 
% dominio de la frecuencia, con la amplitud normalizada y centrada en 0
title('Señal modulada  en el dominio de la frecuencia, f_d = 100[Hz]');
xlabel('frecuencia [Hz]');
ylabel('|F(Hz)|');
axis([9.39*10^5 9.43*10^5 0 0.01])



figure()
subplot(2,1,1);
plot(0:1/fs_up:1,modulada2); %Gráfico de la señal modulante en el 
% dominio de la frecuencia, con la amplitud normalizada y centrada en 0
title('Señal modulada f_d = 1[KHz]');
xlabel('Tiempo [s]');
ylabel('Amplitud [V]');
axis([0 inf -2 2])

subplot(2,1,2);
plot(fshift,fftshift(abs(fft(modulada2)/fs_up))) %Gráfico de la señal modulante en el 
% dominio de la frecuencia, con la amplitud normalizada y centrada en 0
title('Señal modulada  en el dominio de la frecuencia,f_d = 1[kHz]');
xlabel('frecuencia [Hz]');
ylabel('|F(Hz)|');
axis([9.36*10^5 9.47*10^5 0 0.06])



figure()
subplot(2,1,1);
plot(0:1/fs_up:1,modulada3); %Gráfico de la señal modulante en el 
% dominio de la frecuencia, con la amplitud normalizada y centrada en 0
title('Señal modulada f_d = 10[kHz]');
xlabel('Tiempo [s]');
ylabel('Amplitud [V]');
axis([0 inf -2 2])

subplot(2,1,2);
plot(fshift,fftshift(abs(fft(modulada3)/fs_up))) %Gráfico de la señal modulante en el 
% dominio de la frecuencia, con la amplitud normalizada y centrada en 0
title('Señal modulada  en el dominio de la frecuencia, f_d = 10[kHz]');
xlabel('frecuencia [Hz]');
ylabel('|F(Hz)|');
axis([9*10^5 9.82*10^5 0 0.026])



figure()
subplot(2,1,1);
plot(0:1/fs_up:1,modulada4); %Gráfico de la señal modulante en el 
% dominio de la frecuencia, con la amplitud normalizada y centrada en 0
title('Señal modulada f_d = 50[kHz]');
xlabel('Tiempo [s]');
ylabel('Amplitud [V]');
axis([0 inf -2 2])

subplot(2,1,2);
plot(fshift,fftshift(abs(fft(modulada4)/fs_up))) %Gráfico de la señal modulante en el 
% dominio de la frecuencia, con la amplitud normalizada y centrada en 0
title('Señal modulada  en el dominio de la frecuencia, f_d = = 50[kHz]');
xlabel('frecuencia [Hz]');
ylabel('|F(Hz)|');
axis([7.5*10^5 11.32*10^5 0 0.01])


%% Parte 3

ruido1_5 =  ruido(modulada1,5);
ruido1_15 = ruido(modulada1,15);
ruido1_30 = ruido(modulada1,30);

ruido2_5 =  ruido(modulada2,5);
ruido2_15 = ruido(modulada2,15);
ruido2_30 = ruido(modulada2,30);

ruido3_5 =  ruido(modulada3,5);
ruido3_15 = ruido(modulada3,15);
ruido3_30 = ruido(modulada3,30);

ruido4_5 =  ruido(modulada4,5);
ruido4_15 = ruido(modulada4,15);
ruido4_30 = ruido(modulada4,30);


figure()
subplot(4,1,1);
plot(0:1/fs_up:1,modulada1); %Gráfico de la señal modulante en el 
% dominio de la frecuencia, con la amplitud normalizada y centrada en 0
title('Señal modulada f_d = 100[Hz]');
xlabel('Tiempo [s]');
ylabel('Amplitud [V]');
axis([0 inf -2 2])

subplot(4,1,2);
plot(0:1/fs_up:1,ruido1_5); %Gráfico de la señal modulante en el 
% dominio de la frecuencia, con la amplitud normalizada y centrada en 0
title('Señal modulada f_d = 100[Hz] con ruido SNR_{dB}=5');
xlabel('Tiempo [s]');
ylabel('Amplitud [V]');
axis([0 inf -3.5 3.5])

subplot(4,1,3);
plot(0:1/fs_up:1,ruido1_15); %Gráfico de la señal modulante en el 
% dominio de la frecuencia, con la amplitud normalizada y centrada en 0
title('Señal modulada f_d = 100[Hz] con ruido SNR_{dB}=15');
xlabel('Tiempo [s]');
ylabel('Amplitud [V]');
axis([0 inf -3.5 3.5])

subplot(4,1,4);
plot(0:1/fs_up:1,ruido1_30); %Gráfico de la señal modulante en el 
% dominio de la frecuencia, con la amplitud normalizada y centrada en 0
title('Señal modulada f_d = 100[Hz] con ruido SNR_{dB}=30');
xlabel('Tiempo [s]');
ylabel('Amplitud [V]');
axis([0 inf -3.5 3.5])





figure()
subplot(4,1,1);
fshift=(-N/2:N/2-1)*df; %Vector de frecuencias centrado en 0.
plot(fshift,fftshift(abs(fft(modulada1)/fs)))
title('Señal modulada  en el dominio de la frecuencia, f_d = 100[Hz]');
xlabel('frecuencia [Hz]');
ylabel('|F(Hz)|');
axis([9.39*10^5 9.43*10^5 0 0.01])
subplot(4,1,2);
plot(fshift,fftshift(abs(fft(ruido1_5)/fs)))
title('Señal modulada  en el dominio de la frecuencia, f_d = 100[Hz], SNR_{dB}=5');
xlabel('frecuencia [Hz]');
ylabel('|F(Hz)|');
axis([9.39*10^5 9.43*10^5 0 0.01])
subplot(4,1,3);
plot(fshift,fftshift(abs(fft(ruido1_15)/fs)))
title('Señal modulada  en el dominio de la frecuencia, f_d = 100[Hz], SNR_{dB}=15');
xlabel('frecuencia [Hz]');
ylabel('|F(Hz)|');
axis([9.39*10^5 9.43*10^5 0 0.01])
subplot(4,1,4);
plot(fshift,fftshift(abs(fft(ruido1_30)/fs)))
title('Señal modulada  en el dominio de la frecuencia, f_d = 100[Hz], SNR_{dB}=30');
xlabel('frecuencia [Hz]');
ylabel('|F(Hz)|');
axis([9.39*10^5 9.43*10^5 0 0.01])






figure()
subplot(4,1,1);
plot(0:1/fs_up:1,modulada2); %Gráfico de la señal modulante en el 
% dominio de la frecuencia, con la amplitud normalizada y centrada en 0
title('Señal modulada f_d = 1[KHz]');
xlabel('Tiempo [s]');
ylabel('Amplitud [V]');
axis([0 inf -2 2])

subplot(4,1,2);
plot(0:1/fs_up:1,ruido2_5); %Gráfico de la señal modulante en el 
% dominio de la frecuencia, con la amplitud normalizada y centrada en 0
title('Señal modulada f_d = 1[kHz] con ruido SNR_{dB}=5');
xlabel('Tiempo [s]');
ylabel('Amplitud [V]');
axis([0 inf -3.5 3.5])

subplot(4,1,3);
plot(0:1/fs_up:1,ruido2_15); %Gráfico de la señal modulante en el 
% dominio de la frecuencia, con la amplitud normalizada y centrada en 0
title('Señal modulada f_d = 1[kHz] con ruido SNR_{dB}=15');
xlabel('Tiempo [s]');
ylabel('Amplitud [V]');
axis([0 inf -3.5 3.5])

subplot(4,1,4);
plot(0:1/fs_up:1,ruido2_30); %Gráfico de la señal modulante en el 
% dominio de la frecuencia, con la amplitud normalizada y centrada en 0
title('Señal modulada f_d = 1[kHz] con ruido SNR_{dB}=30');
xlabel('Tiempo [s]');
ylabel('Amplitud [V]');
axis([0 inf -3.5 3.5])

figure()
subplot(4,1,1);
fshift=(-N/2:N/2-1)*df; %Vector de frecuencias centrado en 0.
plot(fshift,fftshift(abs(fft(modulada2)/fs)))
title('Señal modulada  en el dominio de la frecuencia, f_d = 1[kHz]');
xlabel('frecuencia [Hz]');
ylabel('|F(Hz)|');
axis([9.36*10^5 9.47*10^5 0 0.06])
subplot(4,1,2);
plot(fshift,fftshift(abs(fft(ruido2_5)/fs)))
title('Señal modulada  en el dominio de la frecuencia, f_d = 1[kHz], SNR_{dB}=5');
xlabel('frecuencia [Hz]');
ylabel('|F(Hz)|');
axis([9.36*10^5 9.47*10^5 0 0.06])
subplot(4,1,3);
plot(fshift,fftshift(abs(fft(ruido2_15)/fs)))
title('Señal modulada  en el dominio de la frecuencia, f_d = 1[kHz], SNR_{dB}=15');
xlabel('frecuencia [Hz]');
ylabel('|F(Hz)|');
axis([9.36*10^5 9.47*10^5 0 0.06])
subplot(4,1,4);
plot(fshift,fftshift(abs(fft(ruido2_30)/fs)))
title('Señal modulada  en el dominio de la frecuencia, f_d = 1[kHz], SNR_{dB}=30');
xlabel('frecuencia [Hz]');
ylabel('|F(Hz)|');
axis([9.36*10^5 9.47*10^5 0 0.06])


figure()
subplot(4,1,1);
plot(0:1/fs_up:1,modulada3); %Gráfico de la señal modulante en el 
% dominio de la frecuencia, con la amplitud normalizada y centrada en 0
title('Señal modulada f_d = 10[KHz]');
xlabel('Tiempo [s]');
ylabel('Amplitud [V]');
axis([0 inf -2 2])

subplot(4,1,2);
plot(0:1/fs_up:1,ruido3_5); %Gráfico de la señal modulante en el 
% dominio de la frecuencia, con la amplitud normalizada y centrada en 0
title('Señal modulada f_d = 10[kHz] con ruido SNR_{dB}=5');
xlabel('Tiempo [s]');
ylabel('Amplitud [V]');
axis([0 inf -3.5 3.5])

subplot(4,1,3);
plot(0:1/fs_up:1,ruido3_15); %Gráfico de la señal modulante en el 
% dominio de la frecuencia, con la amplitud normalizada y centrada en 0
title('Señal modulada f_d = 10[kHz] con ruido SNR_{dB}=15');
xlabel('Tiempo [s]');
ylabel('Amplitud [V]');
axis([0 inf -3.5 3.5])

subplot(4,1,4);
plot(0:1/fs_up:1,ruido3_30); %Gráfico de la señal modulante en el 
% dominio de la frecuencia, con la amplitud normalizada y centrada en 0
title('Señal modulada f_d = 10[kHz] con ruido SNR_{dB}=30');
xlabel('Tiempo [s]');
ylabel('Amplitud [V]');
axis([0 inf -3.5 3.5])

figure()
subplot(4,1,1);
fshift=(-N/2:N/2-1)*df; %Vector de frecuencias centrado en 0.
plot(fshift,fftshift(abs(fft(modulada3)/fs)))
title('Señal modulada  en el dominio de la frecuencia, f_d = 10[kHz]');
xlabel('frecuencia [Hz]');
ylabel('|F(Hz)|');
axis([9*10^5 9.82*10^5 0 0.026])
subplot(4,1,2);
plot(fshift,fftshift(abs(fft(ruido3_5)/fs)))
title('Señal modulada  en el dominio de la frecuencia, f_d = 10[kHz], SNR_{dB}=5');
xlabel('frecuencia [Hz]');
ylabel('|F(Hz)|');
axis([9*10^5 9.82*10^5 0 0.026])
subplot(4,1,3);
plot(fshift,fftshift(abs(fft(ruido3_15)/fs)))
title('Señal modulada  en el dominio de la frecuencia, f_d = 10[kHz], SNR_{dB}=15');
xlabel('frecuencia [Hz]');
ylabel('|F(Hz)|');
axis([9*10^5 9.82*10^5 0 0.026])
subplot(4,1,4);
plot(fshift,fftshift(abs(fft(ruido3_30)/fs)))
title('Señal modulada  en el dominio de la frecuencia, f_d = 10[kHz], SNR_{dB}=30');
xlabel('frecuencia [Hz]');
ylabel('|F(Hz)|');
axis([9*10^5 9.82*10^5 0 0.026])





figure()
subplot(4,1,1);
plot(0:1/fs_up:1,modulada4); %Gráfico de la señal modulante en el 
% dominio de la frecuencia, con la amplitud normalizada y centrada en 0
title('Señal modulada f_d = 50[kHz]');
xlabel('Tiempo [s]');
ylabel('Amplitud [V]');
axis([0 inf -2 2])

subplot(4,1,2);
plot(0:1/fs_up:1,ruido4_5); %Gráfico de la señal modulante en el 
% dominio de la frecuencia, con la amplitud normalizada y centrada en 0
title('Señal modulada f_d = 50[kHz] con ruido SNR_{dB}=5');
xlabel('Tiempo [s]');
ylabel('Amplitud [V]');
axis([0 inf -3.5 3.5])

subplot(4,1,3);
plot(0:1/fs_up:1,ruido4_15); %Gráfico de la señal modulante en el 
% dominio de la frecuencia, con la amplitud normalizada y centrada en 0
title('Señal modulada f_d = 50[kHz] con ruido SNR_{dB}=15');
xlabel('Tiempo [s]');
ylabel('Amplitud [V]');
axis([0 inf -3.5 3.5])

subplot(4,1,4);
plot(0:1/fs_up:1,ruido4_30); %Gráfico de la señal modulante en el 
% dominio de la frecuencia, con la amplitud normalizada y centrada en 0
title('Señal modulada f_d = 50[kHz] con ruido SNR_{dB}=30');
xlabel('Tiempo [s]');
ylabel('Amplitud [V]');
axis([0 inf -3.5 3.5])

figure()
subplot(4,1,1);
fshift=(-N/2:N/2-1)*df; %Vector de frecuencias centrado en 0.
plot(fshift,fftshift(abs(fft(modulada4)/fs)))
title('Señal modulada  en el dominio de la frecuencia, f_d = 50[kHz]');
xlabel('frecuencia [Hz]');
ylabel('|F(Hz)|');
axis([7.5*10^5 11.32*10^5 0 0.01])
subplot(4,1,2);
plot(fshift,fftshift(abs(fft(ruido4_5)/fs)))
title('Señal modulada  en el dominio de la frecuencia, f_d = 50[kHz], SNR_{dB}=5');
xlabel('frecuencia [Hz]');
ylabel('|F(Hz)|');
axis([7.5*10^5 11.32*10^5 0 0.01])
subplot(4,1,3);
plot(fshift,fftshift(abs(fft(ruido4_15)/fs)))
title('Señal modulada  en el dominio de la frecuencia, f_d = 50[kHz], SNR_{dB}=15');
xlabel('frecuencia [Hz]');
ylabel('|F(Hz)|');
axis([7.5*10^5 11.32*10^5 0 0.01])
subplot(4,1,4);
plot(fshift,fftshift(abs(fft(ruido4_30)/fs)))
title('Señal modulada  en el dominio de la frecuencia, f_d = 50[kHz], SNR_{dB}=30');
xlabel('frecuencia [Hz]');
ylabel('|F(Hz)|');
axis([7.5*10^5 11.32*10^5 0 0.01])

%% Parte 4


coseno = cos(2*pi*94.1*10^4*t);
seno   = -1.*sin(2*pi*94.1*10^4*t);

% moduladai ruidoi_dB


%% demodulacion señal 1
modulacion1_1_cos = modulada1 .* coseno;
modulacion1_1_fase = lowpass(modulacion1_1_cos,1000,fs);
modulacion1_1_sen = modulada1 .* seno;
modulacion1_1_cuadratura = lowpass(modulacion1_1_sen,1000,fs);
phi_1_1=atand(modulacion1_1_cuadratura ./ modulacion1_1_fase);
phi_1_1= diff(phi_1_1)./ (2*pi*100);

modulacion1_5_cos = ruido1_5 .* coseno;
modulacion1_5_fase = lowpass(modulacion1_5_cos,1000,fs);
modulacion1_5_sen = ruido1_5 .* seno;
modulacion1_5_cuadratura = lowpass(modulacion1_5_sen,1000,fs);
phi_1_5=atan(modulacion1_5_cuadratura ./ modulacion1_5_fase);
phi_1_5= diff(phi_1_5)./ (2*pi*100);

modulacion1_15_cos = ruido1_15 .* coseno;
modulacion1_15_fase = lowpass(modulacion1_15_cos,1000,fs);
modulacion1_15_sen = ruido1_15 .* seno;
modulacion1_15_cuadratura = lowpass(modulacion1_15_sen,1000,fs);
phi_1_15=atan(modulacion1_15_cuadratura ./ modulacion1_15_fase);
phi_1_15= diff(phi_1_15)./ (2*pi*100);

modulacion1_30_cos = ruido1_30 .* coseno;
modulacion1_30_fase = lowpass(modulacion1_30_cos,200,fs);
modulacion1_30_sen = ruido1_30 .* seno;
modulacion1_30_cuadratura = lowpass(modulacion1_30_sen,200,fs);
phi_1_30=atan(modulacion1_30_cuadratura ./ modulacion1_30_fase);
phi_1_30= diff(phi_1_30)./ (2*pi*100);

figure()
subplot(4,1,1);
plot(0:dt:1-dt,phi_1_1); %Gráfico de la señal modulante en el 
% dominio de la frecuencia, con la amplitud normalizada y centrada en 0
title('Señal demodulada f_d = 100[Hz]');
xlabel('Tiempo [s]');
ylabel('Amplitud [V]');
subplot(4,1,2);
plot(0:dt:1-dt,phi_1_5); %Gráfico de la señal modulante en el 
% dominio de la frecuencia, con la amplitud normalizada y centrada en 0
title('Señal demodulada f_d = 100[Hz], SNR_{dB} = 5');
xlabel('Tiempo [s]');
ylabel('Amplitud [V]');
subplot(4,1,3);
plot(0:dt:1-dt,phi_1_15); %Gráfico de la señal modulante en el 
% dominio de la frecuencia, con la amplitud normalizada y centrada en 0
title('Señal demodulada f_d = 100[Hz], SNR_{dB} = 15');
xlabel('Tiempo [s]');
ylabel('Amplitud [V]');
subplot(4,1,4);
plot(0:dt:1-dt,phi_1_30); %Gráfico de la señal modulante en el 
% dominio de la frecuencia, con la amplitud normalizada y centrada en 0
title('Señal demodulada f_d = 100[Hz], SNR_{dB} = 30');
xlabel('Tiempo [s]');
ylabel('Amplitud [V]');










%% demodulacion señal 2
modulacion2_1_cos = modulada2 .* coseno;
modulacion2_1_fase = lowpass(modulacion2_1_cos,1500,fs);
modulacion2_1_sen = modulada2 .* seno;
modulacion2_1_cuadratura = lowpass(modulacion2_1_sen,1500,fs);
phi_2_1=atan(modulacion2_1_cuadratura ./ modulacion2_1_fase);
phi_2_1= diff(phi_2_1)./ (2*pi*1000);

modulacion2_5_cos = ruido2_5 .* coseno;
modulacion2_5_fase = lowpass(modulacion2_5_cos,1500,fs);
modulacion2_5_sen = ruido2_5 .* seno;
modulacion2_5_cuadratura = lowpass(modulacion2_5_sen,1500,fs);
phi_2_5=atan(modulacion2_5_cuadratura ./ modulacion2_5_fase);
phi_2_5= diff(phi_2_5)./ (2*pi*1000);

modulacion2_15_cos = ruido2_15 .* coseno;
modulacion2_15_fase = lowpass(modulacion2_15_cos,1500,fs);
modulacion2_15_sen = ruido2_15 .* seno;
modulacion2_15_cuadratura = lowpass(modulacion2_15_sen,1500,fs);
phi_2_15=atan(modulacion2_15_cuadratura ./ modulacion2_15_fase);
phi_2_15= diff(phi_2_15)./ (2*pi*1000);

modulacion2_30_cos = ruido2_30 .* coseno;
modulacion2_30_fase = lowpass(modulacion2_30_cos,1500,fs);
modulacion2_30_sen = ruido2_30 .* seno;
modulacion2_30_cuadratura = lowpass(modulacion2_30_sen,1500,fs);
phi_2_30=atan(modulacion2_30_cuadratura ./ modulacion2_30_fase);
phi_2_30= diff(phi_2_30)./ (2*pi*1000);

figure()
subplot(4,1,1);
plot(0:dt:1-dt,phi_2_1); %Gráfico de la señal modulante en el 
% dominio de la frecuencia, con la amplitud normalizada y centrada en 0
title('Señal demodulada f_d = 1[kHz]');
xlabel('Tiempo [s]');
ylabel('Amplitud [V]');
subplot(4,1,2);
plot(0:dt:1-dt,phi_2_5); %Gráfico de la señal modulante en el 
% dominio de la frecuencia, con la amplitud normalizada y centrada en 0
title('Señal demodulada f_d = 1[kHz], SNR_{dB} = 5');
xlabel('Tiempo [s]');
ylabel('Amplitud [V]');
subplot(4,1,3);
plot(0:dt:1-dt,phi_2_15); %Gráfico de la señal modulante en el 
% dominio de la frecuencia, con la amplitud normalizada y centrada en 0
title('Señal demodulada f_d = 1[kHz], SNR_{dB} = 15');
xlabel('Tiempo [s]');
ylabel('Amplitud [V]');
subplot(4,1,4);
plot(0:dt:1-dt,phi_2_30); %Gráfico de la señal modulante en el 
% dominio de la frecuencia, con la amplitud normalizada y centrada en 0
title('Señal demodulada f_d = 1[kHz], SNR_{dB} = 30');
xlabel('Tiempo [s]');
ylabel('Amplitud [V]');




%% demodulacion señal 3
modulacion3_1_cos = modulada3 .* coseno;
modulacion3_1_fase = lowpass(modulacion3_1_cos,12000,fs);
modulacion3_1_sen = modulada3 .* seno;
modulacion3_1_cuadratura = lowpass(modulacion3_1_sen,12000,fs);
phi_3_1=atan(modulacion3_1_cuadratura ./ modulacion3_1_fase);
phi_3_1= diff(phi_3_1)./ (2*pi*10000);

modulacion3_5_cos = ruido3_5 .* coseno;
modulacion3_5_fase = lowpass(modulacion3_5_cos,12000,fs);
modulacion3_5_sen = ruido3_5 .* seno;
modulacion3_5_cuadratura = lowpass(modulacion3_5_sen,12000,fs);
phi_3_5=atan(modulacion3_5_cuadratura ./ modulacion3_5_fase);
phi_3_5= diff(phi_3_5)./ (2*pi*10000);

modulacion3_15_cos = ruido3_15 .* coseno;
modulacion3_15_fase = lowpass(modulacion3_15_cos,12000,fs);
modulacion3_15_sen = ruido3_15 .* seno;
modulacion3_15_cuadratura = lowpass(modulacion3_15_sen,12000,fs);
phi_3_15=atan(modulacion3_15_cuadratura ./ modulacion3_15_fase);
phi_3_15= diff(phi_3_15)./ (2*pi*10000);

modulacion3_30_cos = ruido3_30 .* coseno;
modulacion3_30_fase = lowpass(modulacion3_30_cos,12000,fs);
modulacion3_30_sen = ruido3_30 .* seno;
modulacion3_30_cuadratura = lowpass(modulacion3_30_sen,12000,fs);
phi_3_30=atan(modulacion3_30_cuadratura ./ modulacion3_30_fase);
phi_3_30= diff(phi_3_30)./ (2*pi*10000);

figure()
subplot(4,1,1);
plot(0:dt:1-dt,phi_3_1); %Gráfico de la señal modulante en el 
% dominio de la frecuencia, con la amplitud normalizada y centrada en 0
title('Señal demodulada f_d = 10[kHz]');
xlabel('Tiempo [s]');
ylabel('Amplitud [V]');
subplot(4,1,2);
plot(0:dt:1-dt,phi_3_5); %Gráfico de la señal modulante en el 
% dominio de la frecuencia, con la amplitud normalizada y centrada en 0
title('Señal demodulada f_d = 10[kHz], SNR_{dB} = 5');
xlabel('Tiempo [s]');
ylabel('Amplitud [V]');
subplot(4,1,3);
plot(0:dt:1-dt,phi_3_15); %Gráfico de la señal modulante en el 
% dominio de la frecuencia, con la amplitud normalizada y centrada en 0
title('Señal demodulada f_d = 10[kHz], SNR_{dB} = 15');
xlabel('Tiempo [s]');
ylabel('Amplitud [V]');
subplot(4,1,4);
plot(0:dt:1-dt,phi_3_30); %Gráfico de la señal modulante en el 
% dominio de la frecuencia, con la amplitud normalizada y centrada en 0
title('Señal demodulada f_d = 10[kHz], SNR_{dB} = 30');
xlabel('Tiempo [s]');
ylabel('Amplitud [V]');

%% demodulacion señal 4
modulacion4_1_cos = modulada4 .* coseno;
modulacion4_1_fase = lowpass(modulacion4_1_cos,53000,fs);
modulacion4_1_sen = modulada4 .* seno;
modulacion4_1_cuadratura = lowpass(modulacion4_1_sen,53000,fs);
phi_4_1=atan(modulacion4_1_cuadratura ./ modulacion4_1_fase);
phi_4_1= diff(phi_4_1)./ (2*pi*50000);

modulacion4_5_cos = ruido4_5 .* coseno;
modulacion4_5_fase = lowpass(modulacion4_5_cos,53000,fs);
modulacion4_5_sen = ruido4_5 .* seno;
modulacion4_5_cuadratura = lowpass(modulacion4_5_sen,53000,fs);
phi_4_5=atan(modulacion4_5_cuadratura ./ modulacion4_5_fase);
phi_4_5= diff(phi_4_5)./ (2*pi*50000);

modulacion4_15_cos = ruido4_15 .* coseno;
modulacion4_15_fase = lowpass(modulacion4_15_cos,53000,fs);
modulacion4_15_sen = ruido4_15 .* seno;
modulacion4_15_cuadratura = lowpass(modulacion4_15_sen,53000,fs);
phi_4_15=atan(modulacion4_15_cuadratura ./ modulacion4_15_fase);
phi_4_15= diff(phi_4_15)./ (2*pi*50000);

modulacion4_30_cos = ruido4_30 .* coseno;
modulacion4_30_fase = lowpass(modulacion4_30_cos,53000,fs);
modulacion4_30_sen = ruido4_30 .* seno;
modulacion4_30_cuadratura = lowpass(modulacion4_30_sen,53000,fs);
phi_4_30=atan(modulacion4_30_cuadratura ./ modulacion4_30_fase);
phi_4_30= diff(phi_4_30)./ (2*pi*50000);

figure()
subplot(4,1,1);
plot(0:dt:1-dt,phi_4_1); %Gráfico de la señal modulante en el 
% dominio de la frecuencia, con la amplitud normalizada y centrada en 0
title('Señal demodulada f_d = 50[kHz]');
xlabel('Tiempo [s]');
ylabel('Amplitud [V]');
subplot(4,1,2);
plot(0:dt:1-dt,phi_4_5); %Gráfico de la señal modulante en el 
% dominio de la frecuencia, con la amplitud normalizada y centrada en 0
title('Señal demodulada f_d = 50[kHz], SNR_{dB} = 5');
xlabel('Tiempo [s]');
ylabel('Amplitud [V]');
subplot(4,1,3);
plot(0:dt:1-dt,phi_4_15); %Gráfico de la señal modulante en el 
% dominio de la frecuencia, con la amplitud normalizada y centrada en 0
title('Señal demodulada f_d = 50[kHz], SNR_{dB} = 15');
xlabel('Tiempo [s]');
ylabel('Amplitud [V]');
subplot(4,1,4);
plot(0:dt:1-dt,phi_4_30); %Gráfico de la señal modulante en el 
% dominio de la frecuencia, con la amplitud normalizada y centrada en 0
title('Señal demodulada f_d = 50[kHz], SNR_{dB} = 30');
xlabel('Tiempo [s]');
ylabel('Amplitud [V]');


%% Espectros de frec
figure()
t=0:1/fs:1;
N=length(t); %Largo del vector de Tiempo.
df=1/(N*dt); %Paso de la frecuencia en Hertz.
subplot(4,1,1);
fshift=(-N/2:N/2-2)*df; %Vector de frecuencias centrado en 0.
plot(fshift,fftshift(abs(fft(phi_1_1)/fs)))
title('Señal modulada  en el dominio de la frecuencia, f_d = 100[Hz]');
xlabel('frecuencia [Hz]');
ylabel('|F(Hz)|');
axis([250 3100 0 5*10^-7])
subplot(4,1,2);
plot(fshift,fftshift(abs(fft(phi_1_5)/fs)))
title('Señal modulada  en el dominio de la frecuencia, f_d = 100[Hz], SNR_{dB}=5');
xlabel('frecuencia [Hz]');
ylabel('|F(Hz)|');
axis([250 3100 0 2*10^-8])
subplot(4,1,3);
plot(fshift,fftshift(abs(fft(phi_1_15)/fs)))
title('Señal modulada  en el dominio de la frecuencia, f_d = 100[Hz], SNR_{dB}=15');
xlabel('frecuencia [Hz]');
ylabel('|F(Hz)|');
axis([250 3100 0 2*10^-8])
subplot(4,1,4);
plot(fshift,fftshift(abs(fft(phi_1_30)/fs)))
title('Señal modulada  en el dominio de la frecuencia, f_d = 100[Hz], SNR_{dB}=30');
xlabel('frecuencia [Hz]');
ylabel('|F(Hz)|');
axis([250 3100 0 2*10^-8])


figure()
t=0:1/fs:1;
N=length(t); %Largo del vector de Tiempo.
df=1/(N*dt); %Paso de la frecuencia en Hertz.
subplot(4,1,1);
fshift=(-N/2:N/2-2)*df; %Vector de frecuencias centrado en 0.
plot(fshift,fftshift(abs(fft(phi_2_1)/fs)))
title('Señal modulada  en el dominio de la frecuencia, f_d = 1[kHz]');
xlabel('frecuencia [Hz]');
ylabel('|F(Hz)|');
axis([250 3100 0 5*10^-8])
subplot(4,1,2);
plot(fshift,fftshift(abs(fft(phi_2_5)/fs)))
title('Señal modulada  en el dominio de la frecuencia, f_d = 1[kHz], SNR_{dB}=5');
xlabel('frecuencia [Hz]');
ylabel('|F(Hz)|');
axis([250 3100 0 2*10^-8])
subplot(4,1,3);
plot(fshift,fftshift(abs(fft(phi_2_15)/fs)))
title('Señal modulada  en el dominio de la frecuencia, f_d = 1[kHz], SNR_{dB}=15');
xlabel('frecuencia [Hz]');
ylabel('|F(Hz)|');
axis([250 3100 0 2*10^-8])
subplot(4,1,4);
plot(fshift,fftshift(abs(fft(phi_2_30)/fs)))
title('Señal modulada  en el dominio de la frecuencia, f_d = 1[kHz], SNR_{dB}=30');
xlabel('frecuencia [Hz]');
ylabel('|F(Hz)|');
axis([250 3100 0 2*10^-8])

figure()
t=0:1/fs:1;
N=length(t); %Largo del vector de Tiempo.
df=1/(N*dt); %Paso de la frecuencia en Hertz.
subplot(4,1,1);
fshift=(-N/2:N/2-2)*df; %Vector de frecuencias centrado en 0.
plot(fshift,fftshift(abs(fft(phi_3_1)/fs)))
title('Señal modulada  en el dominio de la frecuencia, f_d = 10[kHz]');
xlabel('frecuencia [Hz]');
ylabel('|F(Hz)|');
axis([250 3100 0 5*10^-9])
subplot(4,1,2);
plot(fshift,fftshift(abs(fft(phi_3_5)/fs)))
title('Señal modulada  en el dominio de la frecuencia, f_d = 10[kHz], SNR_{dB}=5');
xlabel('frecuencia [Hz]');
ylabel('|F(Hz)|');
axis([250 3100 0 2*10^-9])
subplot(4,1,3);
plot(fshift,fftshift(abs(fft(phi_3_15)/fs)))
title('Señal modulada  en el dominio de la frecuencia, f_d = 10[kHz], SNR_{dB}=15');
xlabel('frecuencia [Hz]');
ylabel('|F(Hz)|');
axis([250 3100 0 2*10^-9])
subplot(4,1,4);
plot(fshift,fftshift(abs(fft(phi_3_30)/fs)))
title('Señal modulada  en el dominio de la frecuencia, f_d = 10[kHz], SNR_{dB}=30');
xlabel('frecuencia [Hz]');
ylabel('|F(Hz)|');
axis([250 3100 0 2*10^-9])


figure()
t=0:1/fs:1;
N=length(t); %Largo del vector de Tiempo.
df=1/(N*dt); %Paso de la frecuencia en Hertz.
subplot(4,1,1);
fshift=(-N/2:N/2-2)*df; %Vector de frecuencias centrado en 0.
plot(fshift,fftshift(abs(fft(phi_4_1)/fs)))
title('Señal modulada  en el dominio de la frecuencia, f_d = 50[kHz]');
xlabel('frecuencia [Hz]');
ylabel('|F(Hz)|');
axis([250 3100 0 5*10^-10])
subplot(4,1,2);
plot(fshift,fftshift(abs(fft(phi_4_5)/fs)))
title('Señal modulada  en el dominio de la frecuencia, f_d = 50[kHz], SNR_{dB}=5');
xlabel('frecuencia [Hz]');
ylabel('|F(Hz)|');
axis([250 3100 0 2*10^-10])
subplot(4,1,3);
plot(fshift,fftshift(abs(fft(phi_4_15)/fs)))
title('Señal modulada  en el dominio de la frecuencia, f_d = 50[kHz], SNR_{dB}=15');
xlabel('frecuencia [Hz]');
ylabel('|F(Hz)|');
axis([250 3100 0 2*10^-10])
subplot(4,1,4);
plot(fshift,fftshift(abs(fft(phi_4_30)/fs)))
title('Señal modulada  en el dominio de la frecuencia, f_d = 50[kHz], SNR_{dB}=30');
xlabel('frecuencia [Hz]');
ylabel('|F(Hz)|');
axis([250 3100 0 2*10^-10])