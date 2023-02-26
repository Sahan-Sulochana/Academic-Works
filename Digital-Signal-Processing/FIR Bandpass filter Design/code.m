%%%%%%%%%%%%%General specifications%%%%%%%%%%%%%
A = 2;
B = 3;
C = 7;
A_p = 0.03+(0.01*A); % dB max passband ripple
A_a = 45+B; %dB min stopband attenuation
wp1 = C*100+300; %rad/s lower passband edge
wp2 = C*100+700; %rad/s upper passband edge
wa1 = C*100+150; %rad/s lower stopband edge
wa2 = C*100+800; %rad/s upper stopband edge
ws = 2*((C*100)+1200); %sampling freqency
%%%%%%%%%%%%%Derived specifications%%%%%%%%%%%%%
bt1 = wp1-wa1; %lower transition width
bt2 = wa2-wp2; % upper transisiton width
bt = min(bt1,bt2); %critical transition width
wc1 = wp1-bt/2; % lower cutoff frequency
wc2 = wp2+bt/2; % upper cutoff frequency
T = 2*pi/ws; % sampling period
%%%%%%%%%%%%%Kaiser Window Parameters%%%%%%%%%%%%%
delta_p = (10^(0.05*A_p) - 1)/ (10^(0.05*A_p) + 1); % calculating delta
delta_a = 10^(-0.05*A_a);
delta = min(delta_p,delta_a);
Aa = -20*log10(delta); % Actual stopband attenuation
if Aa<=21 % Calculating alpha
alpha = 0;
elseif Aa>21 && Aa<= 50
alpha = 0.5842*(Aa-21)^0.4 + 0.07886*(Aa-21);
else
alpha = 0.1102*(Aa-8.7);
end
if Aa <= 21 % Calculating D
D = 0.9222;
else
D = (Aa-7.95)/14.36;
end
N = ceil(ws*D/bt +1); % order of the filter
if mod(N,2) == 0
N = N+1;
end
n = -(N-1)/2:1:(N-1)/2; % length of the filter
beta = alpha*sqrt(1-(2*n/(N-1)).^2);

%%%%%%%%%%%%%Generating IoAlpha%%%%%%%%%%%%%
bessellimit = 200;
IoAlpha = 1;
for k = 1:bessellimit
var1 = (1/factorial(k)*(alpha/2).^k).^2;
IoAlpha = IoAlpha + var1;
end
%%%%%%%%%%%%%Generating IoBeta%%%%%%%%%%%%%
IoBeta = 1;
for k = 1:bessellimit
var2 = (1/factorial(k)*(beta/2).^k).^2;
IoBeta = IoBeta + var2;
end
%%%%%%%%%%%%%Obtaining Kaiser Window%%%%%%%%%%%%%
wknt = IoBeta/IoAlpha;
figure
stem(n,wknt)
xlabel('n')
ylabel('Amplitude')
title('Kaiser Window - Time Domain');
%%%%%%%%%%%%%Obtaining Impulse Response%%%%%%%%%%%%%
n_left = -(N-1)/2:-1; %negative part
hnt_left = 1./(n_left*pi).*(sin(wc2*n_left*T)-sin(wc1*n_left*T));
n_right = 1:(N-1)/2; %Positive part
hnt_right = 1./(n_right*pi).*(sin(wc2*n_right*T)-sin(wc1*n_right*T));
hnt0 = 2/ws*(wc2-wc1); %for n=0
hnt = [hnt_left,hnt0,hnt_right];
figure
stem(n,hnt)
xlabel('n')
ylabel('Amplitude')
title(strcat('Filter Response - Rectangular window - Time Domain'));
figure
[h,w] = freqz(hnt);
w = w/T;
h = 20*log10(h);
plot(w,h)
xlabel('Frequency (rad/s)')
ylabel('Magnitude (dB)')
title(strcat('Filter Response - Rectangular Window - Frequency Domain'));

%%%%%%%%%%%%%%%%%%Plotting the Passband(Rectangular)%%%%%%%%%%%%%%%%%%
[h,w] = freqz(hnt);
w = w/T;
h = 20*log10(h);
start = round(length(w)/(ws/2)*wc1);
finish = round((length(w)/(ws/2)*wc2));
wpass1 = w(start:finish);
hpass1 = (h(start:finish));
plot(wpass1,hpass1)
axis([-inf, inf, -1, 1]);
xlabel('Frequency (rad/s)')
ylabel('Magnitude (dB)')
title('Passband - Frequency Domain');
%%%%%%%%%%%%%Applying the Kaiser window to the filter%%%%%%%%%%%%%
filter = hnt.*wknt;
figure
stem(n,filter)
xlabel('n')
ylabel('Amplitude')
title(strcat('Filter Response - Kaiser Window - Time Domain'));
figure
[h,w] = freqz(filter);
w = w/T;
h = 20*log10(abs(h));
plot(w,h)
xlabel('Frequency (rad/s)')
ylabel('Magnitude (dB)')
title(strcat('Filter Response - Kaiser Window - Frequency Domain'));
%%%%%%%%%%%%%Plotting the Passband%%%%%%%%%%%%%
figure
start = round(length(w)/(ws/2)*wc1);
finish = round((length(w)/(ws/2)*wc2));
wpass = w(start:finish);
hpass = (h(start:finish));
plot(wpass,hpass)
axis([-inf, inf, -0.1, 0.1]);
xlabel('Frequency (rad/s)')
ylabel('Magnitude (dB)')
title('Passband - Frequency Domain');

%%%%%%%%%%%%%Input signal generation%%%%%%%%%%%%%
w1 = wc1/2;
w2 = (wc1 + wc2)/2;
w3 = (wc2 + (ws/2))/2;
%Generate discrete signal and evelope
samples = 500;
n1 = 0:1:samples;
n2 = 0:0.1:samples;
x_nT = sin(w1.*n1.*T)+sin(w2.*n1.*T)+sin(w3.*n1.*T);
x_env = sin(w1.*n2.*T)+sin(w2.*n2.*T)+sin(w3.*n2.*T);
%%%%%%%%%%%%%Filtering using frequency domain multiplication%%%%%%%%%%%%%
len_fft = length(x_nT)+length(filter)-1; %length for fft in x dimension
x_fft = fft(x_nT,len_fft);
filter_fft = fft(filter,len_fft);
out_fft = filter_fft.*x_fft; %a shift in time is added here
out = ifft(out_fft,len_fft);
rec_out = out(floor(N/2)+1:length(out)-floor(N/2)); %account for shifting delay
%Ideal Output Signal
ideal_out = sin(w2.*n2.*T);
figure
subplot(2,1,1)
len_fft = 2^nextpow2(numel(n1))-1;
x_fft = fft(x_nT,len_fft);
x_fft_plot = [abs([x_fft(len_fft/2+1:len_fft)]),abs(x_fft(1)),abs(x_fft(2:len_fft/2+1))];
f = ws*linspace(0,1,len_fft)-ws/2;
plot(f,x_fft_plot);
xlabel('Frequency rad/s')
ylabel('Magnitude')
title(strcat(['Input signal',' ','- Frequency Domain']));
%%%%%%%%%%%%%Time domain representation of input signal before filtering%%%%%%
subplot(2,1,2)
stem(n1,x_nT)
axis([0, 50, -inf, inf]);
xlabel('n')
ylabel('Amplitude')
title(strcat(['Input signal',' ','- Time Domain']));
hold on
plot(n2,x_env)
legend('Input signal','Envelope of the Input signal');

%%%%%%%%%%%%%Frequency domain representation of output signal after filtering%
figure
subplot(2,1,1)
len_fft = 2^nextpow2(numel(n1))-1;
xfft_out = fft(rec_out,len_fft);
x_fft_out_plot = [abs([xfft_out(len_fft/2+1:len_fft)]),abs(xfft_out(1)),abs(xfft_out(2:len_fft/2+1))];
f = ws*linspace(0,1,len_fft)-ws/2;
plot(f,x_fft_out_plot); axis tight;
xlabel('Frequency rad/s')
ylabel('Magnitude')
title(strcat(['Output signal',' ','- Frequency Domain']));
%%%%%%%%%%%%%Time domain representation of output signal after filtering%%%%%
subplot(2,1,2)
stem(n1,rec_out)
axis([0, 50, -inf, inf]);
xlabel('n')
ylabel('Amplitude')
title(strcat(['Output signal',' ','- Time Domain']));
hold on
plot(n2,ideal_out)
legend('Output signal','Envelope of the ideal output signal');
fvtool(filter);