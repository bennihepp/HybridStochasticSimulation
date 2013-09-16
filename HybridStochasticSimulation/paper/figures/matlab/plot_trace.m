function p = plot_trace(T, X, indices)
    W = zeros(length(indices), size(X,2));
    i = 1;
    for j = indices
        W(i,:) = X(j,:);
        i = i + 1;
    end
    p = plot(T, W);
    xlabel('time');
    ylabel('copy number');
end
