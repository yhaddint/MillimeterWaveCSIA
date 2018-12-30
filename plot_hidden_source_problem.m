clear;clc;
N = 256;
FFT_size = 256;
freq1 = 0/256*pi;
freq2 = 4/256*pi;
phase1 = exp(1j*0);
phase2 = exp(1j*pi);
sig1 = 10*phase1 * exp(1j*freq1*(0:N-1).');
sig2 = 10*phase2 * exp(1j*freq2*(0:N-1).');
AWGN = (randn(N,1) + 1j*randn(N,1))/sqrt(2)*sqrt(1000);

sig = sig1 + sig2 +AWGN;


%% FFT test
xdata1 = linspace(-pi,pi,FFT_size);
xdata2 = linspace(-pi,pi,FFT_size/8);
ytemp1 = abs(fft(sig(1:FFT_size)));
ytemp2 = abs(fft(sig(1:FFT_size/8)));

figure
subplot(121)
plot(xdata1/pi,20*log10([ytemp1(FFT_size/2:end);ytemp1(1:FFT_size/2-1)]));hold on
grid on
xlabel('Freq. [rad/\pi]')
grid on
ylim([30,80])
xlim([-1,1])
title('Fine FFT')
ylabel('PSD [dB]')

subplot(122)
plot(xdata2/pi,20*log10([ytemp2(FFT_size/16:end);ytemp2(1:FFT_size/16-1)]));hold on
grid on
xlabel('Freq. [rad/\pi]')
grid on
xlim([-1,1])
title('Coarse FFT (8X Wider Window)')
ylim([30,80])
ylabel('PSD [dB]')