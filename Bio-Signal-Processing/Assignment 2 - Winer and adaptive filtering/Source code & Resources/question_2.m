%% MATLAB assignment 2
% Wiener and adaptive filtering
% 
% Sahan Sulochana Hettiarachchi
% 
% 180237G
% 
% BME
%% *2. Adaptive filtering*
% 

clc;
clear all;

% Data Construction
fs = 500; %sampling frequency
N = 5000; % Number of samples
snr = 10; % Signal to Noise Ratio

t = linspace(0, 5, N); % Time vector
s = sawtooth(2*pi*2*t(1 : N), 0.5); % Sawtooth signal
n1 = 0.2*sin(2*pi*50*t(1 : N/2)); % Sinusoidal noise with 50 Hz
n2 = 0.3*sin(2* pi*100*t(N/2 + 1 : N));% Sinusoidal noise with 100 Hz
nwg = s - awgn(s, snr, 'measured'); % 10 dB Gaussian white noise

noise_added_signal = s + nwg + [n1 n2]; % Noise added signal  

figure('Name', 'Signals using for Adaptive Filtering')
subplot(2,1,1)
plot(t, s)
title('Sawtooth wave (ideal signal)')
subplot(2,1,2)
plot(t, noise_added_signal)
title('Noise added signal (noisy signal)')
linkaxes()
%%
% arbitary constants
a = 1.611; 
pi_1 = pi*(1/6); 
pi_2 = pi*(1/2);

% Non-stationary noise 
n_50r = 0.2*sin(2*pi*50*t(1 : N/2) + pi_1); 
n_100r = 0.3*sin(2* pi*100*t(N/2 + 1 : N) + pi_2); 

% Constructed signal
r_n = a*(nwg + [n_50r n_100r]);
%% 2.1. LMS method
%%
%% LMS_method Algorithm

mu_1 = 0.006112; %Converging factor
M = 12; % Filter order

% Appliying LMS function 
[error, ~, ~] = LMS_method(noise_added_signal, r_n, mu_1, M);

% Plotting tht figures
figure;
subplot(4,1,1)
plot(t,s)
title('Desired Signal sawtooth_sig(n)')
subplot(4,1,2)
plot(t,noise_added_signal)
title('Noise Corrupted Signal')
subplot(4,1,3)
plot(t,error)
title(['Filtered Signal using LMS_method M=10 u=0.006161'])
subplot(4,1,4)
plot(t,abs(error-s))
title('Absolute Error')

LMS_error = error;

% Calculating Mean Squared Error
M_range = 15;
MSE = NaN(M_range, 100);

max_lambda = 20*M_range*((noise_added_signal*noise_added_signal')/ length(noise_added_signal));
mu = linspace(0.001, 2/ max_lambda, 100);

for M = 1:M_range
    for i = 1:100
        [error, ~, ~] = LMS_method(noise_added_signal, r_n, mu(i), M);
        MSE(M,i) = immse(error,s);
    end
end
%% 
% 

% Plotting figures
M = 1:M_range;
surf(mu, M, MSE)
title('MSE Variation when Adaptive Filtering using LMS') 
colorbar
xlabel('mu'), ylabel('M - Order'),zlabel('MSE');
colormap('jet');
% Find LMS with minimum MSE

[ms, ls] = min(MSE,[],2);
[min_MSE, min_m] = min(ms);
min_lambda = ls(min_m)*(2/ max_lambda - 0.001)/100 + 0.001;
disp(['Minimum Error = ' num2str(min_MSE) ' at M = ' num2str(min_m) ' and mu = ' num2str(min_lambda)]);
%% 2.2. RLS method
%%
% fastest RLS_method filter

lambda = 0.996; % 0 < lamda <= 1
M = 15; % Order

[error_2, ~, ~] = RLS_method(noise_added_signal, r_n, lambda, M);

figure;
subplot(4,1,1)
plot(t, s)
title('Desired Signal sawtooth_sig(n)')
subplot(4,1,2)
plot(t, noise_added_signal)
title('Noise Corrupted Signal')
subplot(4,1,3)
plot(t, error_2);
title('Filtered Signal using RLS_method M=15 lambda = 0.996');
subplot(4,1,4)
plot(t, abs(error_2 - s'))

RLS_error = error_2;

% Calculating MSE 
M_range = 15;
MSE = NaN(M_range,100);
lambda = linspace(0.9,1,100);

for M=1:M_range
    for i = 1:100
        error = RLS_method(noise_added_signal, r_n, lambda(i), M);
        MSE(M,i) = immse(error', s);
    end
end
%% 
% 
%%
% Plotting figures 
figure;
surf(lambda,(1:M_range), MSE)
colorbar
title('MSE Variation when Adaptive Filtering using RLS')
xlabel('lambda'), ylabel('M - Order'), zlabel('MSE');
colormap('jet')

%% Find RLS with min MSE
[ms,ls] = min(MSE,[],2);
[min_MSE,m_min] = min(ms);
min_lambda = ls(m_min)*(0.01)/100 + 0.9;
disp(['Minimum Error = ' num2str(min_MSE) ' at M = ' num2str(m_min) ' and lambda = ' num2str(min_lambda)])
%%
% Comparing LMS_method and RLS_method

figure
subplot(2,1,1)
plot(t, abs(LMS_error - s));
title(['Error convergence using the LMS algorithm \mu = ' num2str(mu_1) ' M = 15' ]);
xlabel('Time(sawtooth_sig)');
ylabel('Voltage (mV)');
grid on
subplot(2,1,2)
plot(t, abs(RLS_error' - s))
title(['Filtered Signal of the ANC filter using the RLS algorithm \lambda = ' num2str(lambda) ' M = 15']),
xlabel('Time(sawtooth_sig)');
ylabel('Voltage (mV)');
grid on
%%
% Adaptive FIltering ECG_sig signal

%loading the ECG signal
load('idealECG.mat')
ECG_signal = idealECG - mean(idealECG);

fs = 500; % Sampling frequency
N = length(ECG_signal); % Length of the signal
t = linspace(0,N/fs,N);                  

signal = ECG_signal;                                 
n1 = 0.2*sin(2*pi*50*t(1 : N/2));                     
n2 = 0.3*sin(2* pi*100*t(N/2+1 : N));                

% Adding White Gaussian Noise with SNR=10
nwg = signal - awgn(signal, 10,'measured');  

% Noise added signal
noise_added_signal = signal + nwg + [n1 n2];         
%%
% Signal generation

% arbitary constants
a = 1.611; 
pi_1 = pi*(1/6); 
pi_2 = pi*(1/2);

% Non-stationary noise
n_50r = 0.2*sin(2*pi*50*t(1 : N/2) + pi_1); 
n_100r = 0.3*sin(2*pi*100*t(N/2 + 1 : N) + pi_2); 

% Non-stationary noise added signal
r_n = a*(nwg + [n_50r n_100r]);
%% 
% *LMS method applying for ECG signal*
%%
mu = 0.006112;
M = 12;
[error, ~, ~] = LMS_method(noise_added_signal, r_n, mu, M);
%% 
% 
% 
% *RLS method applying for ECG signal*
%%
lambda = 0.996; % 0 < lamda <= 1
M = 12;
[error_2, y2, w2] = RLS_method(noise_added_signal, r_n, lambda, M);
%% 
% *Plotting graphs*
%%
figure;
subplot(6,1,1)
plot(t, signal)
title('Desired Signal sawtooth signal(n)')
subplot(6,1,2)
plot(t, noise_added_signal)
title('Noise Corrupted Signal')
subplot(6,1,3)
plot(t, error)
title(['Filtered Signal using LMS method M=12 u=0.006112'])
subplot(6,1,4)
plot(t, abs(error - signal))
title('Absolute Error | LMS')
subplot(6,1,5)
plot(t,error_2);
title(['Filtered Signal using RLS method M=12 lambda = 0.996']);
subplot(6,1,6)
plot(t,abs(error_2 - signal'))
title('Absolute Error | RLS')
linkaxes()