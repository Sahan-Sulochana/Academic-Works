function y_hat = weiner_td(input_signal, w_matrix)

    order = length(w_matrix);
    y_hat = input_signal;
    %applying the filter for the signal with the weight matrix
    for i = 1: length(input_signal) - order
        y_hat(i) = input_signal(i : i + order - 1) * w_matrix;
    end
end