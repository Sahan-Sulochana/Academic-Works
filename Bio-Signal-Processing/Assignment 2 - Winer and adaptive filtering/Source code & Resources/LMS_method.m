function [error_vector, new_signal, weight_matrix] = LMS_method(noise_sig, r_n, mu, M)

    Nr = length(r_n);
    Nx = length(noise_sig);
    
    % Initialize outputs and required matrices
    
    Rr = zeros(M,1);% signal delay matrix
    weight_vector = zeros(M,1);% adaptive weights vector
    weight_matrix = zeros(Nx,M);% adaptive weights matrix
    new_signal = zeros(1, Nx);
    error_vector = zeros(1, Nx);% Error vector

    if (Nr <= M)  
        disp('error: Signal length is less than the filter order');
        return; 
    end
    if (Nr ~= Nx)  
        disp('error: Input signal and reference signal are different in length?');
        return; 
    end

    lambda_max = 20*M*((noise_sig*noise_sig')/length(noise_sig));
    
    if (mu > 2/lambda_max) 
        disp(['mu is too large' num2str(mu) ' /' num2str(lambda_max)]);
        return
    end


    for k = 1:Nx
        Rr(1) = r_n(k);% sig_R(n)
        
        % new_sig(n) = weight_mat(n)'.sig_R(n) new_sig gives the error
        new_signal(k) = weight_vector'*Rr;
        
        % err_vect(n) = noisy_sig(n) - weight_mat(n)'.sig_R(n) = noisy_sig(n) - new_sig(n)
        error_vector(k) = noise_sig(k) - new_signal(k);
        
        % weight_mat(n+1) = weight_mat(n) + 2*miu*err_vect(n)*sig_R(n) - Widrow Hoff LMS algorthm
        weight_vector = weight_vector + 2*mu*error_vector(k)*Rr; 
        
        % Store the weights vector in the weights matrix
        weight_matrix(k,:) = weight_vector; 
        
        % Delay back the reference signal window by one sample
        Rr(2:M) = Rr(1:M-1); 
    end

end