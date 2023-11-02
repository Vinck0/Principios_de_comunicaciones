function [y,n]= ruido(x,SNR_dB)
    %[y,n]=ruido(x,SNR) agrega un vector de ruido AWGN  a señal ’x’ 
    % y es la señal resultante y n el vector de ruido
    L=length(x);
    SNR = 10^(SNR_dB/10); 
    S=sum(abs(x).^2)/(L); 
    N0=S/SNR; %Densidad espectral del ruido
    noiseSigma = sqrt(N0);%Standard deviation for AWGN Noise 
    n = noiseSigma*randn(1,L);%computed noise
    y = x + n;
end