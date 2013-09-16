function plot_distribution(W, nbins)
    if nargin < 3
        nbins = 10;
    end
    hist(W, nbins);
    xlabel('copy numbers');
    ylabel('count');
end
