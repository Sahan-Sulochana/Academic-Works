% This function is used to reconstruct the signal using waveltes coefficent for level10
% Wavelet : name of the wavelet (Harr or DB9)

% A : Apporiximaion coefficent array
% A : detailed coefficents array at each level
function reconstructed = wave_reconstruction(wavelet, A,D1,D2,D3,D4,D5,D6,D7,D8,D9,D10)

    %calculting previous wavelts
    A_9 = idwt(A,D10,wavelet);
    A_8 = idwt(A_9,D9,wavelet);
    A_7 = idwt(A_8,D8,wavelet);
    A_6 = idwt(A_7,D7,wavelet);
    A_5 = idwt(A_6,D6,wavelet);
    A_4 = idwt(A_5,D5,wavelet);
    
     if wavelet == "haar"
       A_3 = idwt(A_4,D4,wavelet); 
     else
        A_3 = idwt(A_4(1:79),D4,wavelet); 
     end
    A_2 = idwt(A_3,D3,wavelet);
    A_1 = idwt(A_2,D2,wavelet);
    reconstructed = idwt(A_1,D1,wavelet);
    
end