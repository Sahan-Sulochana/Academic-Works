clearvars; close all; clc;

fs = 250;           % Sample frequency
N = 3000;           % Data length
t = (-N:N)/fs;      % Time scale
s = 0.01:0.1:2;     % Values of scaling factor

% Maxican hat wavelet
%...


% Generating spectra of wavelets
Fwavelt = fft(wavelt)/length(wavelt);
hz = linspace(0,fs/2,floor(length(wavelt)/2)+1);
plot(hz,2*abs(Fwavelt(1:length(hz))))
xlabel('Frequency (Hz)'), ylabel('Amplitude')

% % Ploting the spectrogram
% h = pcolor();
% set(h, 'EdgeColor', 'none');
% colormap jet
% xlabel('Time (s)')
% ylabel('Scale')
% title('Spectrogram')