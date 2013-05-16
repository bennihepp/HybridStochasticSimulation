function [t, x] = simulate_SSA(x0, t0, t1, prop_func, stoch_matrix)

    % iteration variable i
    i = 1;

    % set initial time and state
    t = t0;
    x = x0;

    while t <= t1
        % draw sample for transcription and degradation reaction
        %t = exprnd(1 / (k + X(i) * y));
        q = prop_func(x);
        qsum = sum(q);
        tau = exprnd(1 / qsum);
        u = rand(1);
        t = t + tau;
        w = 0;
        for l = 1:length(q)
            w = w + q(l) / qsum;
            if u < w
                x = x + stoch_matrix(l,:)';
                break;
            end
        end
        i = i + 1;
    end

end
