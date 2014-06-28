function [alpha, beta, gamma] ...
    = search_multiscale_alpha_beta_gamma(N, stoch_matrix, x0, prop_rates)

    r0 = size(stoch_matrix, 1);
    s0 = size(stoch_matrix, 2);

    javaaddpath('JaCoP-3.2.jar');
    javaaddpath('bin');
    cp_species_balance = CPSpeciesBalance();
    cp_species_balance.setStochMatrix(stoch_matrix);

    min_alpha = max(round(log(x0) / log(N)) - 1, 0);
    max_alpha = max(round(log(x0) / log(N)), 0) + 1;
    min_beta = round(log(prop_rates) / log(N)) - 1;
    max_beta = round(log(prop_rates) / log(N)) + 1;
    gamma = 0;

    %cp_species_balance.set_structure(mu, beta_plus, eta, beta_minus);
    x = cp_species_balance.search_solution(min_alpha, max_alpha, ...
        min_beta, max_beta, gamma);

    alpha = x(1:s0);
    beta = x(s0+1:end);

%     alpha = zeros(s0,1);
%     beta = zeros(r0,1);
%     alpha = [1; 1; 0; 0];
%     beta = [-1; -1];
%     gamma = 0;
%     holds_ts = check_time_scale_constraints( ...
%                 alpha, beta, gamma, omega, beta_set);
%     holds_be = check_balance_equation( ...
%                 alpha, beta, gamma, mu, beta_plus, eta, beta_minus);
%     holds_ts
%     holds_be
%     all(holds_ts)
%     all(holds_be)

end

function [holds] ...
    = check_time_scale_constraints(alpha, beta, gamma, omega, beta_set)

    s0 = length(alpha);
    r0 = length(beta);

    holds = zeros(s0,1);
    for s = 1:s0
        max_term = -inf;
        for k = 1:size(omega{s},1)
            term = beta(beta_set{s}(k)) + omega{s}(k,:) * alpha;
            if term > max_term
                max_term = term
            end
        end
        holds(s) = gamma <= alpha(s) - max_term;
    end

end

function [holds] = check_balance_equation( ...
    alpha, beta, gamma, mu, beta_plus, eta, beta_minus)

    s0 = length(alpha);
    r0 = length(beta);

    holds = zeros(s0,1);
    for s = 1:s0
        max_plus_term = -inf;
        for k = 1:size(mu{s},1)
            term = beta(beta_plus{s}(k)) + mu{s}(k,:) * alpha;
            if term > max_plus_term
                max_plus_term = term
            end
        end
        max_minus_term = -inf;
        for k = 1:size(eta{s},1)
            term = beta(beta_minus{s}(k)) + eta{s}(k,:) * alpha;
            if term > max_minus_term
                max_minus_term = term
            end
        end
        holds(s) = (max_plus_term == max_minus_term) ...
                    && max_plus_term ~= -inf;
    end

end
