Comparison of Analytic Wavelet Transform and Analytic Signal Coefficients
This example shows how the analytic wavelet transform of a real signal approximates the corresponding analytic signal. 
This is demonstrated using a sine wave. If you obtain the wavelet transform of a sine wave using an analytic wavelet and extract the wavelet coefficients at a scale corresponding to the frequency of the sine wave, the coefficients approximate the analytic signal. For a sine wave, the analytic signal is a complex exponential of the same frequency.
Create a sinusoid with a frequency of 50 Hz.
t = 0:.001:1;
x = cos(2*pi*50*t);

Obtain its continuous wavelet transform using an analytic Morse wavelet and the analytic signal. You must have the Signal Processing Toolbox™ to use hilbert.
[wt,f] = cwt(x,1000,'voices',32,'ExtendSignal',false);
analytsig = hilbert(x);

Obtain the wavelet coefficients at the scale closest to the sine wave's frequency of 50 Hz.
[~,idx] = min(abs(f-50));
morsecoefx = wt(idx,:);

Compare the real and imaginary parts of the analytic signal with the wavelet coefficients at the signal frequency.
figure;
plot(t,[real(morsecoefx)' real(analytsig)']);
title('Real Parts'); 
ylim([-2 2]); grid on;
legend('Wavelet Coefficients','Analytic Signal','Location','SouthEast');
xlabel('Time'); ylabel('Amplitude');

figure;
plot(t,[imag(morsecoefx)' imag(analytsig)']);
title('Imaginary Parts'); 
ylim([-2 2]); grid on;
legend('Wavelet Coefficients','Analytic Signal','Location','SouthEast');
xlabel('Time'); ylabel('Amplitude');
cwt uses L1 normalization and scales the wavelet bandpass filters to have a peak magnitude of 2. The factor of 1/2 in the above equation is canceled by the peak magnitude value.
The wavelet transform represents a frequency-localized filtering of the signal. Accordingly, the CWT coefficients are less sensitive to noise than are the Hilbert transform coefficients. 

Add highpass noise to the signal and reexamine the wavelet coefficients and the analytic signal.
y = x + filter(1,[1 0.9],0.1*randn(size(x)));
analytsig = hilbert(y);
[wt,f] = cwt(y,1000,'voices',32,'ExtendSignal',0);
morsecoefy = wt(idx,:);

figure;
plot(t,[real(analytsig)' x']);
legend('Analytic Signal','Original Signal');
grid on;
xlabel('Time'); ylabel('Amplitude');
ylim([-2 2])

figure;
plot(t,[real(morsecoefy)' x']);
legend('Wavelet Coefficients','Original Signal');
grid on;
xlabel('Time'); ylabel('Amplitude');
ylim([-2 2])
Copyright 2012 The MathWorks, Inc.
