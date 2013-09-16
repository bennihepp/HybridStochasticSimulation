function p = plot_mean(T, X, indices)
    Q = squeeze(mean(X, 1));
    p = plot_trace(T, Q, indices);
end
