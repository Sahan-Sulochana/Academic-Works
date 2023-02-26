function [y_hat, Wf] = weiner_fd(temp, noise, in_signal)

    N = length(in_signal);
    
    % Power of template signal
    Syy = abs(fft(temp,N*2-1)).^2;
    
    % Power of noise signal
    SNN = abs(fft(noise,N*2-1)).^2;

    % Fourier transform of the signal
    SXf  = fft(in_signal,N*2-1);
    
    Wf = Syy./(Syy + SNN); % Weinier filter in frequency domain
    S_Yhat = Wf.*SXf;      % Estimated signal in frequency domain
    
    y_hat_td = (ifft(S_Yhat));% Time domain signal
    y_hat = y_hat_td(1: length(in_signal));% Since FFT being two sided

end