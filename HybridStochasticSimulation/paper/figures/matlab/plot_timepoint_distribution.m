function t = plot_timepoint_distribution(t, T, X, indices, labels, nbins)
    k = find(T >= t, 1);
    t = T(k);
    Q = squeeze(X(:,:,k));
    W = zeros(size(Q,1), length(indices));
    i = 1;
    for j = indices
        W(:,i) = Q(:,j);
        i = i + 1;
    end
    i = 1;
    cmap = colormap('Lines');
    for j = indices
        subplot(length(indices), 1, i);
        W = Q(:,j);
        if nargin < 7
            plot_distribution(W);
        else
            plot_distribution(W, nbins);
        end
        c = cmap(i,:);
        h = findobj(gca, 'Type', 'patch');
        set(h, 'FaceColor', c)
        legend(labels(i));
        i = i + 1;
    end
end
