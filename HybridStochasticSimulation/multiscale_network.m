% Note that this function is limited to unary and binary reactions
function [z, ms_stoch_matrix, ...
          ms_produce_stoch_matrix, ms_consume_stoch_matrix, ...
          ms_prop_func, ms_evo_func, ms_info] = ...
    multiscale_network(x, produce_stoch_matrix, ...
                       consume_stoch_matrix, prop_rates, ...
                        N, gamma, alpha, beta)

    stoch_matrix = produce_stoch_matrix - consume_stoch_matrix;

    r0 = size(produce_stoch_matrix,1);
    s0 = size(produce_stoch_matrix,2);

    % Compute scaled state
	z = N.^(-alpha) .* x;
    % Compute scaled propensity rates
    ms_prop_rates = N.^(-beta) .* prop_rates;

    % Compute rho for each reaction,
    % build expressions for the combinatorial choices
    % for each reaction in stochastic and deterministic form and
    % build an expression for the scaled propensities
    rho = gamma+beta;
    ms_prop_helper_func_str = '@(z) [';
    ms_prop_choices_str = {};
    ms_evo_choices_str = {};
    for r = 1:r0
        % q will contain all the consumed species of this reaction
        q = [];
        for s = 1:s0
            nu = consume_stoch_matrix(r,s);
            if nu > 0
                q(end+1) = s;
                switch nu
                    case 1
                        rho(r) = rho(r) + alpha(s);
                    case 2
                        rho(r) = rho(r) + 2*alpha(s);
                    otherwise
                        error(['only unary or binary ' ...
                                'reactions are allowed!']);
                end
            end
        end
        % Here we build the combinatorial choices depending on the
        % type of reaction
        switch length(q)
            case 0
                % Constitutive reaction
                ms_prop_choices_str{r} = '1';
                ms_evo_choices_str{r} = ms_prop_choices_str{r};
            case 1
                if consume_stoch_matrix(r, q(1)) == 2
                    % Binary reaction of the same species
                    ms_prop_choices_str{r} = ['1/2*z(', num2str(q(1)), ...
                        ')*(z(', num2str(q(1)), ...
                        ') - N^-alpha(', num2str(q(1)), '))'];
                    ms_evo_choices_str{r} = ['1/2*z(', num2str(q(1)), ...
                        ')^2'];
                else
                    % Unary reaction
                    ms_prop_choices_str{r} = ['z(', ...
                        num2str(q(1)), ')'];
                    ms_evo_choices_str{r} = ms_prop_choices_str{r};
                end
            case 2
                % Binary reaction
                ms_prop_choices_str{r} = ['z(', ...
                    num2str(q(1)), ')*z(', num2str(q(2)), ')'];
                ms_evo_choices_str{r} = ms_prop_choices_str{r};
            otherwise
                error(['only unary or binary ' ...
                        'reactions are allowed!']);
        end
        ms_prop_helper_func_str = [ms_prop_helper_func_str, ...
            '', ms_prop_choices_str{r}, '; '];
    end
    ms_prop_helper_func_str = [ms_prop_helper_func_str, ']'];

    % Compute the scaling of the individual reaction terms.
    % Outside scaling is the scaling outside of the Poisson process
    outside_scaling = repmat(-alpha', [r0,1]);
    neg_outside_scaling = -outside_scaling;
    % Inside scaling is the scaling inside of the Poisson process
    inside_scaling = repmat(rho, [1,s0]);

    % Here we evaluate the scaled propensity expression to a
    % Matlab function
    ms_prop_helper_func = eval(ms_prop_helper_func_str);
    ms_prop_func = @(z) (N.^rho) .* ms_prop_rates ...
                        .* ms_prop_helper_func(z);

    % Compute a scaled stochiometry matrix
    ms_produce_stoch_matrix = produce_stoch_matrix .* (N.^outside_scaling);
    ms_consume_stoch_matrix = consume_stoch_matrix .* (N.^outside_scaling);
    ms_stoch_matrix = stoch_matrix .* (N.^outside_scaling);

    % Build an expression for the deterministic evolution function
    % in matrix form (one expression per species and reaction)
    ms_evo_helper_func_str = '@(t, z) [';
    for r = 1:r0
        str = '';
        for s = 1:s0
            if stoch_matrix(r,s) ~= 0
                str = [str, num2str(stoch_matrix(r,s)), ...
                        '*ms_prop_rates(', num2str(r), ')*', ...
                        ms_evo_choices_str{r}, ', '];
            else
                str = [str, '0', ', '];
            end
        end
        ms_evo_helper_func_str = [ms_evo_helper_func_str, str, '; '];
    end
    ms_evo_helper_func_str = [ms_evo_helper_func_str, ']'];
    % Here we evaluate the built expression to a Matlab function
    ms_evo_matrix_func = eval(ms_evo_helper_func_str);

    % TODO: Differentiate between producing and consuming reactions?
    % Now we check the convergence of the scaled reaction terms.
    % We compare outside and inside scaling to determine the convergence.
    % Convergence to zero
    die_out_mask = (neg_outside_scaling > inside_scaling) ...
                    & (stoch_matrix ~= 0);
    % Convergence to deterministic term
    deterministic_mask = (neg_outside_scaling == inside_scaling) ...
                            & (neg_outside_scaling > 0) ...
                            & (stoch_matrix ~= 0);
    % Convergence to stochastic term
    stochastic_mask = (neg_outside_scaling == inside_scaling) ...
                        & (neg_outside_scaling == 0) ...
                        & (stoch_matrix ~= 0);
    % No convergence (i.e. exploding terms)
    exploding_mask = (neg_outside_scaling < inside_scaling) ...
                        & (stoch_matrix ~= 0);

    % Now we mask the scaled propensity function, the deterministic
    % evolution function and the scaled stochiometry matrix depending on
    % the convergence of each term.
    % Here we only keep propensities and stochiometries that have a
    % corresponding stochastic term
    ms_prop_func_mask = any(stochastic_mask, 2);
    masked_prop_func = @(z) ms_prop_func(z) .* ms_prop_func_mask;
    masked_produce_stoch_matrix = zeros(size(ms_produce_stoch_matrix));
    masked_produce_stoch_matrix(stochastic_mask) ...
        = ms_produce_stoch_matrix(stochastic_mask);
    masked_consume_stoch_matrix = zeros(size(ms_consume_stoch_matrix));
    masked_consume_stoch_matrix(stochastic_mask) ...
        = ms_consume_stoch_matrix(stochastic_mask);
    masked_stoch_matrix = zeros(size(ms_stoch_matrix));
    masked_stoch_matrix(stochastic_mask) = ms_stoch_matrix(stochastic_mask);
    % Here we only keep evolution terms that have a corresponding
    % deterministic term
    masked_evo_func = @(t, z) sum(ms_evo_matrix_func(t, z) ...
                                .* deterministic_mask, 1)';

    % We return the masked variables/functions
    ms_stoch_matrix = masked_stoch_matrix;
    ms_produce_stoch_matrix = masked_produce_stoch_matrix;
    ms_consume_stoch_matrix = masked_consume_stoch_matrix;
    ms_prop_func = masked_prop_func;
    ms_evo_func = masked_evo_func;
    ms_info = struct;
    ms_info.terms = struct;
    ms_info.terms.die_out = sum(die_out_mask(:));
    ms_info.terms.deterministic = sum(deterministic_mask(:));
    ms_info.terms.stochastic = sum(stochastic_mask(:));
    ms_info.terms.exploding = sum(exploding_mask(:));
    assert(ms_info.terms.die_out + ms_info.terms.deterministic ...
            + ms_info.terms.stochastic + ms_info.terms.exploding ...
            == sum(sum(stoch_matrix ~= 0)));
    ms_info.masks = struct;
    ms_info.masks.die_out = die_out_mask;
    ms_info.masks.deterministic = deterministic_mask;
    ms_info.masks.stochastic = stochastic_mask;
    ms_info.masks.exploding = exploding_mask;
    ms_info.scaled_propensity_rates = ms_prop_rates;
    ms_info.rho = rho;
    ms_info.inside_scaling = inside_scaling;
    ms_info.outside_scaling = outside_scaling;

end
