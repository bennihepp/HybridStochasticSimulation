function [T, R, X] = simulate_PDMP(x0, t0, t1, prop_func, stoch_matrix, ...
                                evo_func)

    N = 1000;
    dt = (t1 - t0) / 1000.0;

    % T will hold time points, X will hold state values
    % and R will hold reaction indices
    T = zeros(1,N);
    R = zeros(1,N);
    X = zeros(N,length(x0));

    % set initial time and state
    t = t0;
    x = x0;
    x = vertcat(x, [0; 0]);
    ext_evo_func = @(t, x) vertcat(evo_func(t, x(1:end-2)), ...
                                    [sum(prop_func(x(1:end-2))); 0]);
    ode_options = odeset('Events', @ode_events);

    T(1) = t0;
    R(1) = 0;
    X(1,:) = x0;

    % index variable i
    i = 1;
    % iteration variable
    j = 1;

    while t < t1
        x(end - 1) = 0;
        x(end) = -log(rand(1));
        % evolve deterministic part until first discrete reaction
        %tspan = [t, t1];
        tspan = t:dt:t1;
        if length(tspan) < 2
            tspan = [t, t1];
        end
        ode_sol = ode45(ext_evo_func, tspan, x, ode_options);
        % draw sample for transcription and degradation reaction
        %t = exprnd(1 / (k + X(i) * y));
        if i + length(ode_sol.x) + 1 > length(T)
            T = horzcat(T, zeros(1,N));
            R = horzcat(R, zeros(1,N));
            X = vertcat(X, zeros(N,length(x0)));
        end

        T(i+1:i+length(ode_sol.x)-2) = ode_sol.x(2:end-1);
        R(i+1:i+length(ode_sol.x)-2) = 0;
        X(i+1:i+length(ode_sol.x)-2,:) = ode_sol.y(1:end-2,2:end-1)';
        i = i + length(ode_sol.x)-2;
%         T = horzcat(T, ode_sol.x(2:end-1));
%         X = vertcat(X, ode_sol.y(1:end-2,2:end-1)');
        t = ode_sol.x(end);
        x = ode_sol.y(:,end);
        q = prop_func(x(1:end-2));
        qsum = sum(q);
        u = rand(1);
        w = 0;
        r = 0;
        for l = 1:length(q)
            w = w + q(l) / qsum;
            if u < w
                r = l;
                x(1:end-2) = x(1:end-2) + stoch_matrix(r,:)';
                break;
            end
        end
        T(i+1) = t;
        R(i+1) = r;
        X(i+1,:) = x(1:end-2)';
        i = i + 1;
        j = j + 1;
    end

    T = T(1:i);
    R = R(1:i);
    X = X(1:i,:);

end

function [value, isterminal, direction] = ode_events(t, x)
    value = x(end) - x(end - 1);
    isterminal = 1;
    direction = 0;
end
