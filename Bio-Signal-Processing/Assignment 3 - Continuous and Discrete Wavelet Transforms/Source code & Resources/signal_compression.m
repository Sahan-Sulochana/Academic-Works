function signal_compression(x, num_level, wavelet_type, thresh_percent, sig_name)

[c, l] = wavedec(x, num_level, wavelet_type);

%% Coefficients sortion in decending order
c_sorted = sort(abs(c(:)),'descend');

% Plotting the sorted coefficients
figure('Name',['Sorted ' wavelet_type ' Wavelet Coefficients of ' sig_name ' - in descending order'])
stem(c_sorted);
xlim([0, length(c_sorted)])
title(['Sorted ' wavelet_type ' Wavelet Coefficients of ' sig_name ' - in descending order']);

%%
cumulative_energy = 0;
num_of_selected_coef = 0;

% Total Energy
total_energy = sum((c_sorted).^2);

for j=1:length(c_sorted)
    cumulative_energy = cumulative_energy + (c_sorted(j)).^2;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if (round(cumulative_energy/total_energy,2) == thresh_percent/100)
        num_of_selected_coef = j;
        break;
    end
end

%% Printing number of coefficients according to the requirement
disp(['Number coefficients required to represent 99% of the energy of the signal = ' num2str(num_of_selected_coef)]);

compression_ratio = length(x)/num_of_selected_coef;
disp(['Compression Ratio = ' num2str(compression_ratio)]);

threshold = c_sorted(num_of_selected_coef);


c_selected = c;
%% Removing noise components by assigning zeros
for k = 1:length(c_selected)
    if (abs(c_selected(k)) < threshold)
        c_selected(k) = 0;
    end
end

% Reconstruct the signal with the remaining coefficients
x_reconst = waverec(c_selected, l, wavelet_type);

% Plotting the reconstructed signal
len_y = length(x);
N = 1:1:len_y;

figure ('Name', [sig_name ' reconstructed with ' wavelet_type]);
plot(N, x_reconst);
xlim([0 len_y]);
title([sig_name ' reconstructed with ' wavelet_type]), xlabel('Samples(n)'), ylabel('Amplitude');

%% Calculating the RMSE between the original signal and the reconstructed signal
rmse = sqrt(immse(x, x_reconst));

%Printing the RMSE value
disp(['RMSE of ' sig_name ' reconstructed with ' wavelet_type ' wavelet = ' num2str(rmse)]);

% Original and reconstructed signal comparison
figure('Name',['Comparing original and reconstructed ' sig_name ' with ' wavelet_type ])

plot(N, x, 'r', N, x_reconst, 'b')
xlim([0 len_y]);
title(['Comparing original and reconstructed ' sig_name ' with ' wavelet_type ]), xlabel('Samples(n)'), ylabel('Amplitude');
legend(sig_name, ['Reconstructed ' sig_name])
end