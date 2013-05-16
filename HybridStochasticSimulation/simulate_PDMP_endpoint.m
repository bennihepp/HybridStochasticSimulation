function [t, x] = simulate_PDMP_endpoint(x0, t0, t1, prop_func, stoch_matrix, ...
                                evo_func, kgoff, kgon, kr, yr, kp, yp)

    % iteration variable i
    i = 1;

    t = t0;
    x = x0;
    x = vertcat(x, [0; 0]);
    ext_evo_func = @(t, x) vertcat(evo_func(t, x(1:end-2)), ...
                                    [sum(prop_func(x(1:end-2))); 0]);
    %ext_evo_func = @(t, x) [0; 0; 0; kp*x(3) - yp*x(4); ...
    %                       kgoff*x(2) + kgon*x(1) + kr*x(2) + yr*x(3); ...
    %                       0];
    ode_options = odeset('Events', @ode_events);

    while t < t1
        x(end - 1) = 0;
        x(end) = -log(rand(1));
        % evolve deterministic part until first discrete reaction
        ode_sol = ode45(ext_evo_func, [t, t1], x, ode_options);
        % draw sample for transcription and degradation reaction
        %t = exprnd(1 / (k + X(i) * y));
        t = ode_sol.x(end);
        x = ode_sol.y(:,end);
        q = prop_func(x(1:end-2));
        qsum = sum(q);
        u = rand(1);
        w = 0;
        for l = 1:length(q)
            w = w + q(l) / qsum;
            if u < w
                break;
            end
        end
        x(1:end-2) = x(1:end-2) + stoch_matrix(l,:)';
        i = i + 1;
    end
    x = x(1:end-2);

end

function [value, isterminal, direction] = ode_events(t, x)
    value = x(end) - x(end - 1);
    isterminal = 1;
    direction = 0;
end
