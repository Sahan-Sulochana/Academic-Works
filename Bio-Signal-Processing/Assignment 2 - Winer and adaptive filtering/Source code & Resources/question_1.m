%% MATLAB assignment 2
% Wiener and adaptive filtering
% 
% Sahan Sulochana Hettiarachchi
% 
% 180237G
% 
% BME
% 
% 
%% 1. Weiner filtering (on stationary signals)
% 1.1. Discrete-time domain implementation of the Wiener filter

load('idealECG.mat'); %loading the ECG file
L = length(idealECG); %length of the ECG signal
N = 1:L;

fs = 500; %sampling frequency
t = linspace(0,length(idealECG)*(1/fs),L);

%ideal ECG signal
ideal_ECG_sig = idealECG - mean(idealECG);

nwg = awgn(ideal_ECG_sig,10,'measured'); %Additive white gaussian Noise
time  = linspace(0, L-1, L)*(1/fs);
n50 = 0.2*sin(2*pi*50*time);
noiseECG = nwg + n50;

%Input signal
x = noiseECG;
%% 
% *Original ECG signal*
%%
figure;
plot(t,ideal_ECG_sig);
xlabel('Time(s)');
ylabel('Magnitude(mV)');
title('The original ECG Signal');
%% 
% *Zoomed ECG signal*
%%
figure;
plot(t,ideal_ECG_sig);
xlabel('Time(s)');
ylabel('Magnitude(mV)');
xlim(1.5:2.5);
%values are obtained by manually
line([1.876 1.876],[-1 2],'Color','red','LineStyle','--');%one of T wave peak
line([2.043 2.043],[-1 2],'Color','red','LineStyle','--');%next T wave peak
%84 samples
title('The original ECG Signal');
%% 
% *Noise ECG signal*
%%
% extracting the noise from the T wave to P wave
figure;
plot(t,x);
xlabel('Time(s)');
ylabel('Magnitude(mV)');
xlim(1.5:2.5);
%values are obtained by manually
line([2.032 2.032],ylim,'Color','red','LineStyle','--');%T wave peak
line([2.072 2.072],ylim,'Color','red','LineStyle','--');%P wave peak
title('The Input Signal');
%% 
% *Single beat of ideal ECG signal*
%%
y1_start = floor(1.876*fs); %starting index of the selecting signal array
y1_stop = floor(2.043*fs); %ending index of the selecting signal array
y1 = ideal_ECG_sig(y1_start:y1_stop); %Selected ideal signal
t_y1 = 0:(1/fs):(y1_stop-y1_start)*(1/fs);

figure;
plot(t_y1,y1);
xlabel('Time(s)');
ylabel('Magnitude(mV)');
%% 
% *Iso electric signal *
%%
x1_start = floor(2.032*fs); %starting index of the selecting signal array
x1_stop = floor(2.072*fs); %ending index of the selecting signal array
x1 = [x(x1_start:x1_stop) x(x1_start:x1_stop) x(x1_start:x1_stop) x(x1_start:x1_stop)]; %Selected signal
xt1 = 0:length(x1)-1;
figure;
plot(xt1,x1);
xlabel('Time(s)');
ylabel('Magnitude(mV)');
title('Isoelectric segment of x(n)');
%% 
% Let take an arbitary filter order as 15
%%
order = 10; % arbitary filter order
noise_signal = noiseECG-ideal_ECG_sig; %noise signal
weight_mat = weight_vector(y1, x1, order);
disp(weight_mat');
y_hat = weiner_td(noiseECG, weight_mat);
figure;
plot(t, ideal_ECG_sig, t, noiseECG, t, y_hat)
legend('Ideal ECG','Noisy ECG','Weiner Filtered ECG (order = 10)')
xlabel('Time (s)'); 
ylabel('Voltage (mV)');
title('Weiner Filtering of noise added ECG (M = 10)');
%% 
% *The optimum filter order and its coefficients*
%%
%Selecting the range as 50 for finding optimised M(filter order)
range = 50;

%Mean Square Error vector
MSE = NaN(1,range);

for M = 2: range
    weight_mat = weight_vector(y1, x1, M);
    y_hat = weiner_td(noiseECG(938:1021),weight_mat); 
    MSE(M) = immse(y_hat, ideal_ECG_sig(938:1021));
end

figure;
plot(MSE)
hold on 
[min_MSE,min_M] = min(MSE); %minimum MSE is marked
scatter(min_M, min_MSE)
title('Weiner Filters MSE vs Different Filter Orders')
xlabel('Filter Order');
ylabel('MSE');

% Optimum Order
optimum_M = min_M;

%Optimum weight vector
optimum_M_weight_mat = weight_vector(y1, x1, optimum_M);

%Optimum filter applied signal
optimum_M_y_hat = weiner_td(noiseECG, optimum_M_weight_mat);

figure;
plot(t, idealECG, t, noiseECG, t, optimum_M_y_hat)
legend('Ideal ECG','Noisy ECG','Weiner Filtered ECG (Optimum order)')
title('Weiner Filtering noise corrupted ECG with Optimum order')
xlabel('Time (s)');
ylabel('Voltage (mV)');
%%
fvtool(optimum_M_weight_mat);
% Plotting the spectrum of yi(n), n(n), x(n) and y_hat(n)

%PSD of noisy ECG signal
[psd_noiseEEG, f1_noisyEEG] = periodogram(noiseECG, [], [], fs);

%PSD of noise signal
[psd_noise,f2_noise] = periodogram(noise_signal, [], [], fs);

%PSD of ideal ECG signal
[psd_ideal_ECG,f3_ideal_ECG] = periodogram(ideal_ECG_sig, [], [], fs);

%PSD of optimal filtered signal
[psd_y_hat,f4_y_hat] = periodogram(optimum_M_y_hat, [], [], fs);

figure('Name','PSD')
semilogy(f1_noisyEEG, psd_noiseEEG, f2_noise, psd_noise, f3_ideal_ECG, psd_ideal_ECG, f4_y_hat, psd_y_hat);
legend('Noise corrupted signal','Noise','Required Signal','Optimum Wiener Filtered Signal')
title('Power Spectral Density'),xlabel('Frequency (Hz)'),ylabel('dB/Hz');
%% 
% *Part 2*
% 
% *Linear model of the ECG signal*
%%
y2 = linearModel(y1);
figure;
plot([0:83],y2);
xlabel('Time(s)');
ylabel('Magnitude(mV)');
title('Linear model of the ideal single beat ECG');
%% 
% 
%%
order = 10; % arbitary filter order
noise_signal = noiseECG-ideal_ECG_sig; %noise signal
weight_mat = weight_vector(y2, x1, order);
disp(weight_mat');%obtain the weight matrix
y_hat = weiner_td(noiseECG, weight_mat);
figure;
plot(t, ideal_ECG_sig, t, noiseECG, t, y_hat)
legend('Ideal ECG','Noisy ECG','Weiner Filtered ECG using linear model(order = 10)')
xlabel('Time (s)'); 
ylabel('Voltage (mV)');
title('Weiner Filtering noise corrupted ECG M = 10');
%% 
% **
% 
% *Optimum order*
%%
%Selecting the range as 50 for finding optimised M(filter order)
range = 50;

%Mean Square Error vector
MSE = NaN(1,range);

for M = 2: range
    weight_mat = weight_vector(y2, x1, M);
    y_hat = weiner_td(noiseECG(938:1021),weight_mat); 
    MSE(M) = immse(y_hat, ideal_ECG_sig(938:1021));
end

figure;
plot(MSE)
hold on 
[min_MSE,min_M] = min(MSE); %minimum MSE is marked
scatter(min_M, min_MSE)
title('Weiner Filters MSE vs Different Filter Orders(For linear model)')
xlabel('Filter Order');
ylabel('MSE');

% Optimum Order
optimum_M = min_M;

%Optimum weight vector
optimum_M_weight_mat = weight_vector(y2, x1, optimum_M);

%Optimum filter applied signal
optimum_M_y_hat = weiner_td(noiseECG, optimum_M_weight_mat);

figure;
plot(t, idealECG, t, noiseECG, t, optimum_M_y_hat)
legend('Ideal ECG','Noisy ECG','Weiner Filtered ECG (Optimum order)')
title('Weiner Filtering noise corrupted ECG with Optimum order(For linear model)')
xlabel('Time (s)');
ylabel('Voltage (mV)');
%%
fvtool(optimum_M_weight_mat);
% Plotting the spectrum of yi(n), n(n), x(n) and y_hat(n)

%PSD of noisy ECG signal
[psd_noiseEEG, f1_noisyEEG] = periodogram(noiseECG, [], [], fs);

%PSD of noise signal
[psd_noise,f2_noise] = periodogram(noise_signal, [], [], fs);

%PSD of ideal ECG signal
[psd_linear_ECG,f3_ideal_ECG] = periodogram(ideal_ECG_sig, [], [], fs);

%PSD of optimal filtered signal
[psd_y_hat,f4_y_hat] = periodogram(optimum_M_y_hat, [], [], fs);

figure('Name','PSD')
semilogy(f1_noisyEEG, psd_noiseEEG, f2_noise, psd_noise, f3_ideal_ECG, psd_linear_ECG, f4_y_hat, psd_y_hat);
legend('Noise added signal','Noise','Ideal ECG Signal','Optimum Wiener Filtered Signal(For linear model)')
title('Power Spectral Density');
xlabel('Frequency (Hz)');
ylabel('dB/Hz');
%% 1.2. Frequency domain implementation of the Wiener filter
%%
[y_hat_ft, W_f] = weiner_fd(ideal_ECG_sig, noise_signal, noiseECG);
figure;
plot(time, ideal_ECG_sig, time, noiseECG, time, y_hat_ft);
legend('Ideal ECG','Noisy ECG','Weiner Filtered ECG');
title('Weiner Filtering Frequency domain');
xlabel('Time (s)');
ylabel('Voltage (mV)');
%%
%Plotting the PSD

%PSD of noisy ECG signal
[psd_noisy_sig, f1_noisy_sig] = periodogram(noiseECG, [], [], fs);

%PSD of noise signal
[psd_noise,f2_noise] = periodogram(noise_signal, [], [], fs);

%PSD of ideal ECG signal
[psd_ideal,f3_ideal] = periodogram(ideal_ECG_sig, [], [], fs);

%PSD of frequency domain filtered signal
[psd_yhat_ft, f4_yhat_ft] = periodogram(y_hat_ft, [], [], fs);

figure('Name','PSD Weiner Filter Frequency domain')
semilogy(f1_noisy_sig, psd_noisy_sig, f2_noise, psd_noise, f3_ideal, psd_ideal, f4_yhat_ft, psd_yhat_ft);
legend('Noise added signal','Noise','Ideal ECG Signal','Frequency domain Wiener Filtered Signal')
title('Power Spectral Density');
xlabel('Frequency (Hz)');
ylabel('dB/Hz');
%%
% Comparing Frequency domain and Time domain weiner filter

figure('Name', 'Comparing Frequency domain and Time domain weiner filter');
plot(time, ideal_ECG_sig, 'r', time, optimum_M_y_hat,'g', time, y_hat_ft, 'b')
legend('Ideal ECG','Filtered by Optimum Time Domain derived Weiner ', 'Filtered by Freq Domain derived Weiner Filter')
title('Comparing Frequency domain and Time domain weiner filter')
xlabel('Time (s)'), ylabel('Voltage (mV)');

%Mean Square Error calculation for time domain derivation
MSE_td = immse(optimum_M_y_hat, ideal_ECG_sig);

%Mean Square Error calculation for frequency domain derivation
MSE_fd = immse(y_hat_ft, ideal_ECG_sig);

disp('Mean Square Error (Time domain)');
disp(MSE_td);
disp('Mean Square Error (Frequency domain)');
disp(MSE_fd);
%% *1.3. Effect on non-stationary noise on the Wiener filtering*
%%
%non-stationary noise Generation

f1 = 50;
f2 = 100;

t_p1 = time(1 : floor(length(time)/2));
t_p2 = time(floor(length(time)/ 2) + 1 : end);

n50_p1 = 0.2*sin(2*pi*f1*t_p1);
n100_p1 = 0.3*sin(2*pi*f2*t_p2);

non_stationary_noise = [n50_p1 n100_p1];
non_stationary_noise_sig = nwg + non_stationary_noise;
%% 
% filtering the non stationary noise added signal with the derived frequency 
% domain weiner filter

N = length(non_stationary_noise_sig);  
S_Xf  = fft(non_stationary_noise_sig, N*2-1);

% Signal estimate from observation and using Wiener filter
S_Yhat = W_f.* S_Xf; 

% Time domain conversion
y_hat_time = ifft(S_Yhat);              
y_hat_non_statationary = y_hat_time(1 : N);


figure;
plot(time, ideal_ECG_sig, time, y_hat_non_statationary, time, y_hat_ft)
xlim([time(1),time(end)])
legend('Ideal ECG','Non-Stationary Noise - Filtered','Stationary Noise - Filtered')
title('Non-Stationary Noise Comparison after filtering with Weiner Freq')
xlabel('Time (s)');
ylabel('voltage(mV)');
%%
noise_signal_2 = non_stationary_noise_sig - ideal_ECG_sig;

%Plotting PSD

%PSD of non-stationary noise added signal
[psd_noisy_sig, f1_noisy_sig] = periodogram(non_stationary_noise_sig, [], [], fs);

%PSD of non-stationary noise signal 
[psd_noise,f2_noise] = periodogram(noise_signal_2, [], [], fs);

%PSD of original ECG signal
[psd_ideal,f3_ideal] = periodogram(ideal_ECG_sig, [], [], fs);

%PSD of weiner filtered (Frequency domain) signal
[psd_yhat_fd_ns, f4_yhat_freq_ns] = periodogram(y_hat_non_statationary, [], [], fs);

figure('Name','PSD Weiner Filter Frequency domain | With Non Sattionary Noise')
semilogy(f1_noisy_sig, psd_noisy_sig, f2_noise, psd_noise, f3_ideal, psd_ideal, f4_yhat_freq_ns, psd_yhat_fd_ns);
legend('Non stationary Noise added signal','Non Sattionary Noise','Ideal ECG Signal','Frequency domain Wiener Filtered Signal')
title('Power Spectral Density | With Non Sattionary Noise');
xlabel('Frequency (Hz)');
ylabel('dB/Hz');
%% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
%