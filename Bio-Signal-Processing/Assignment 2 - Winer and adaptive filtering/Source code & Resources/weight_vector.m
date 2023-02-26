function W = weight_vector(y, est_noise, order)

    yyT = 0;%autocorelation 
    nnT = 0;%autocorelation 
    Yy = 0;%crosscorrelation
    y_vect = zeros(order,1);
    n_vect = zeros(order,1);
    L = length(y);%length of the signal

    for i=1:L
        
        y_vect(1) = y(i); 
        n_vect(1) = est_noise(i);

        yyT = yyT + toeplitz(autocorr(y_vect, order-1));
        nnT = nnT + toeplitz(autocorr(n_vect, order-1));

        Yy = Yy + y_vect*y(i);

        % shifting the delay 
        y_vect(2:order) = y_vect(1 : order-1);
        n_vect(2:order) = n_vect(1 : order-1);
    end

    yyT = yyT.*mean(y.^2);
    nnT = nnT.*mean(est_noise.^2);

    autocorr_Y = yyT./ (L - order); %Ryy
    autocorr_N = nnT./ (L - order); %Rxx
    theta_Yy = Yy./ (L-order); %RYy
    
    autocorr_X = autocorr_Y + autocorr_N; %Rxx = Ryy + Rnn
    W = autocorr_X\theta_Yy;
    %W = inv(autocorr_X)*theta_Yy; %W = inverse(Rxx)* RYy

end