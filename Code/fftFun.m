function [frequency, magnitude, phase] = fftFun(x, fs)

% Fast Fourier Transform
y = fft(x);

m = length(y);
frequencyAxis = (0:m)*fs/m; 
fftAmplitude = abs(y);
fftPhase = angle(y);

frequency = [-(frequencyAxis(m)-frequencyAxis(m/2:m)) frequencyAxis(1:m/2)];
magnitude = [fftAmplitude(m/2:m) fftAmplitude(1:m/2)];
phase = [fftPhase(m/2:m) fftPhase(1:m/2)];


