% Copyright (C) 2021 Stella
% Author: Styliani Karagianni
% Created: 2021-02-15
% Final Project Arrythmia Detection; DBSP 2020-2021

clear all;

% Load signal
val=load('100m.mat');
ECGSignal=struct2array(val);
Fs=360;
t=(0:length(ECGSignal)-1)/Fs;
plot(t,ECGSignal);  %initial signal

% Peaks detection
figure;
hold on;
h = [-1 1]; % dervative filter
d = conv(ECGSignal,h);
maxima = [];
for z=1:length(ECGSignal)
    if ( (d(z) <= 0) & (d(z+1) > 0) ) 
        maxima(z) = ECGSignal(z);
    else
        maxima(z) = 0;
    end
end
plot(maxima,'ro');
plot(ECGSignal);
axis([0 10000 -1050 1050]);  
hold off;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Noise filtering

x=ECGSignal; 
N=length(x);    % number of samples
n=1:N;          % sample index
t=n*(1/Fs);     % t of samples

% Extra power line noise (if wanted)
% pn = 10*cos(2*pi*50*t);
% x=x+pn; 
% figure;
% plot(t,x);

% Spectrum visualization before noise filtering (if wanted)
% figure;
% [frequency,magnitude,phase]=fftFun(x,Fs);
% plot(frequency,magnitude);

X=fft(x);

%breathing noise
cutoffHz = 0.5;
cutoffk = int16(cutoffHz*N/Fs);
X(1:cutoffk)=0;
X(N-cutoffk:N)=0;

%electricity noise
cutoff1Hz = 49;
cutoff1k = int16(cutoff1Hz*N/Fs);
cutoff2Hz = 51;
cutoff2k = int16(cutoff2Hz*N/Fs);
X(cutoff1k:cutoff2k)=0;
X(N-cutoff2k:N-cutoff1k)=0;

%muscle noise
cutoff3Hz = 15;
cutoff3k = int16(cutoff3Hz*N/Fs);
X(1:cutoff3k)=0;
X(N-cutoff3k:N)=0;

% Return in time field after process
y=ifft(X);

% Spectrum visualization of the new signal
figure;
[frequency,magnitude,phase]=fftFun(y,Fs);
plot(frequency,magnitude);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Heart Rate
distances = [];
countZeros = 0;
countPeaks = 0;
for z=1:length(maxima)
    if (maxima(z)>0)
        countPeaks = countPeaks + 1;
        distances(countPeaks) = countZeros;
        countZeros = 0;
    else
        countZeros = countZeros + 1;
    end
end

Ts = 1/Fs;
times = distances * Ts;
bpms = 60 * 1./times

figure;
plot(bpms(2:length(bpms)));
