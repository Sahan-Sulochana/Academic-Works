
function [error, new_signal, weight_matrix] = RLS_method(noise_signal, R, lambda, L)

    %Noise Signal length    
    N_signal = length(noise_signal);
    
    N_r = length(R);

    % Initialize outputs

    I = eye(L);
    alpha = 0.01;
    p = alpha * I;

    % signal delay matrix
    xx = zeros(L,1);
    
    % adaptive weights vector
    weight_vect = zeros(L,1);
    
    % adaptive weights matrix
    weight_matrix = zeros(L,N_signal); 
    
    new_signal = zeros(N_signal,1);
    
    % Error vector
    error = zeros(N_signal,1); 
    
    if (N_r <= L)  
        print('error: Signal length is less than the filter order');
        return; 
    end
    if (N_r ~= N_signal)  
        print('error: Input signal and reference signal are different in length?');
        return; 
    end

    for n = 1:N_signal
        % R(n)
        xx(1) = R(n);                               
        k = (p * xx) ./ (lambda + xx' * p * xx);
        
        % new_sig(n) = weight_mat(n)'.R(n)
        new_signal(n) = xx'*weight_vect;
        
        % e(n) = noisy_sig(n) - weight_mat(n)'.R(n) = noisy_sig(n) - new_sig(n)
        error(n) = noise_signal(n) - new_signal(n);
        
        % weight_mat(n+1) = weight_mat(n) + k**e(n)
        weight_vect = weight_vect + k * error(n);     
        p = (p - k * xx' * p) ./ lambda;
        
        % Store the weights vector in the weights matrix
        weight_matrix(:,n) = weight_vect;
        
        % Delay back the reference signal window by one sample
        xx(2:L) = xx(1:L-1);                        
    end

end