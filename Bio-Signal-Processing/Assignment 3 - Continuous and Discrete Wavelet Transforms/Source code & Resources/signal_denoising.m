function signal_denoising(x_n, y_n, num_levels, wavelet_name, thresh_val, sig_name)

[c, l] = wavedec(y_n, num_levels, wavelet_name);

%% Coefficients sortion in decending order
c_sorted = sort(abs(c(:)),'descend');

% Plotting the sorted coefficients
figure('Name',['Sorted ' wavelet_name ' Wavelet Coefficients of ' sig_name ' - in descending order'])
stem(c_sorted);
xlim([0, length(c_sorted)])
title(['Sorted ' wavelet_name ' Wavelet Coefficients of ' sig_name ' - in descending order']);


c_selected = c;

%% Removing noise components by assigning zeros
for k = 1:length(c_selected)
    if (abs(c_selected(k)) < thresh_val)
        c_selected(k) = 0;
    end
end

% Reconstruct the signal with the remaining coefficients
x_reconst = waverec(c_selected, l, wavelet_name);

%% Plotting the reconstructed signal 
len_y = length(y_n);
N = 1:1:len_y;

figure ('Name', [sig_name ' reconstructed with ' wavelet_name]);
plot(N, x_reconst);
xlim([0 len_y]);
title([sig_name ' reconstructed with ' wavelet_name]), xlabel('Samples(n)'), ylabel('Amplitude');

%% Calculating the RMSE between the original signal and the reconstructed signal
error = x_n - x_reconst;   
rmse = sqrt(sum(abs(error).^2)/length(error));

%Printing the RMSE value
disp(['RMSE of ' sig_name ' reconstructed with ' wavelet_name ' wavelet = ' num2str(rmse)]);

%% Original and reconstructed signal comparison
figure('Name',['Comparing original and reconstructed ' sig_name ' with ' wavelet_name ])
plot(N, x_n, 'r', N, x_reconst, 'b')
xlim([0 len_y]);
title(['Original and reconstructed signal comparison ' sig_name ' with ' wavelet_name ]);
xlabel('Samples(n)');
ylabel('Amplitude');
legend(sig_name, ['Reconstructed ' sig_name])
end