function [T, R, X] = simulate_SSA(x0, t0, t1, prop_func, stoch_matrix)

    N = 1000;

    % T will hold time points, X will hold state values
    % and R will hold reaction indices
    T = zeros(1,N);
    R = zeros(1,N);
    X = zeros(N,length(x0));

    % iteration variable i
    i = 1;

    % set initial time and state
    t = t0;
    x = x0;

    T(1) = t0;
    R(1) = 0;
    X(1,:) = x0;

    while t < t1
        % draw sample for transcription and degradation reaction
        %t = exprnd(1 / (k + X(i) * y));
        q = prop_func(x);
        qsum = sum(q);
        tau = exprnd(1 / qsum);
        u = rand(1);
        t = t + tau;
        w = 0;
        r = 0;
        for l = 1:length(q)
            w = w + q(l) / qsum;
            if u < w
                r = l;
                x = x + stoch_matrix(r,:)';
                break;
            end
        end
        if i + 1 > length(T)
            T = horzcat(T, zeros(1,N));
            R = horzcat(R, zeros(1,N));
            X = vertcat(X, zeros(N,length(x0)));
        end
        T(i+1) = t;
        R(i+1) = r;
        X(i+1,:) = x';
        i = i + 1;
    end

    T = T(1:i);
    R = R(1:i);
    X = X(1:i,:);

end
