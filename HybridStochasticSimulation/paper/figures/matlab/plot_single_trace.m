function p = plot_single_trace(T, X, n, indices)
    Q = squeeze(X(n,:,:));
    plot_trace(T, Q, indices);
end
