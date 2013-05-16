function example_hybrid_network_simulation()

    %rng(101);

%     simple_crystallization_example_network;
%     gamma = 0;
%     alpha = [1; 1; 0; 0];
%     beta = [-1; -1];

%     two_species_birth_death_example_network;
%     gamma = 0;
%     alpha = [0; 0];
%     beta = [0; 0; 0; 0];

    regulated_transcription_example_network;
    gamma = 0;
    alpha = [1; 1; 0; 0; 0; 0];
    beta = [-1; -2; -1; -1; -1; 0; -3; -2; -1; 0];

    javaaddpath('JaCoP-3.2.jar');
    javaaddpath('jgrapht-jdk1.6.jar');
    javaaddpath('commons-math3-3.2.jar');
    javaaddpath('bin');
%     g = ReactionNetworkGraph();
%     g.setStochMatrix(stoch_matrix);
%     stronglyConnCompLists = g.stronglyConnectedComponentsAsLists();
%     stronglyConnComp = {};
%     for i = 0:stronglyConnCompLists.size()-1
%         stronglyConnComp{i+1} = [];
%         compList = stronglyConnCompLists.get(i);
%         for j = 0:compList.size()-1
%             stronglyConnComp{i+1}(end+1) = compList.get(j) + 1;
%         end
%     end

%     [alpha, beta, gamma] = search_multiscale_alpha_beta_gamma( ...
%         N, stoch_matrix, x0, prop_rates)

    [z0, ms_stoch_matrix, ~, ~, ms_prop_func, ms_evo_func, ms_info] ...
        = multiscale_network(x0, production_stoch_matrix, ...
                             consumption_stoch_matrix, prop_rates, ...
                             N, gamma, alpha, beta);

    % initial time t0 and end time t1
    t0 = 0;
    t1 = 100;

    if ms_info.terms.exploding > 0
        error('Some terms in the scaled process are exploding');
    end
    [num2str(ms_info.terms.die_out), ' terms in the scaled process ' ...
        'are dying out']
    [num2str(ms_info.terms.deterministic), ' terms in the scaled ' ...
        'process are deterministic']
    [num2str(ms_info.terms.stochastic), ' terms in the scaled ' ...
        'process are stochastic']

    stoch_matrix_pdmp = ms_stoch_matrix;
    prop_func_pdmp = ms_prop_func;
    evo_func_pdmp = ms_evo_func;

    simulateDistributions = 0;
    simulateTrajectories = 1;
    simulateTrajectoryDistributions = 0;
    runSSA = 1;
    runPDMP = 1;
    runWithJava = 0;

    if simulateDistributions
        M = 3;

        % simulate distributions
        'SSA'
        Xssa = zeros(M,s0);
        for i = 1:M
            i
%             [t, x] = simulate_SSA_endpoint(x0, t0, t1, prop_func, ...
%                                     stoch_matrix);
%             Xssa(i,:) = x;
            [T, R, X] = simulate_SSA(x0, t0, t1, prop_func, stoch_matrix);
            if T(end) > t1
                Xssa(i,:) = X(end-1,:);
            else
                Xssa(i,:) = X(end,:);
            end
        end
        'PDMP'
        Xpdmp = zeros(M,s0);
        for i = 1:M
            i
%             [t, x] = simulate_PDMP_endpoint(x0, t0, t1, prop_func_pdmp, ...
%                                     stoch_matrix, evo_func_pdmp);
%             Xpdmp(i,:) = x;
            [T, R, X] = simulate_PDMP(x0, t0, t1, prop_func_pdmp, ...
                                    stoch_matrix_pdmp, evo_func_pdmp);
            if T(end) > t1
                Xpdmp(i,:) = X(end-1,:);
            else
                Xpdmp(i,:) = X(end,:);
            end
        end

        figure(1);
        subplot(2, 1, 1);
        %plot_trajectories(Tssa, Xssa, species_str, 'SSA Trajectories');
        plot_distribution(t1, Xssa, species_str, plot_scaling, ...
                            'SSA distributions');
        subplot(2, 1, 2);
        %plot_trajectories(Tpdmp, Xpdmp, species_str, 'PDMP Trajectories');
        plot_distribution(t1, Xpdmp, species_str, plot_scaling, ...
                            'PDMP distributions');
    end

    if simulateTrajectories
        % simulate trajectories
        if runSSA
            'SSA'
            [Tssa, Rssa, Xssa] = simulate_SSA(x0, t0, t1, prop_func, ...
                                    stoch_matrix);
            ['SSA: simulated ', num2str(sum(Rssa > 0)), ' reactions']
        end
        if runPDMP
            'PDMP'
            if runWithJava
                [Tpdmp, Zpdmp, reactionTimes, reactions] ...
                    = simulate_PDMP_java( ...
                        x0, t0, t1, ...
                        production_stoch_matrix, ...
                        consumption_stoch_matrix, ...
                        prop_rates, ...
                        N, gamma, alpha, beta, ...
                        1000);
                Xpdmp = Zpdmp;
                ['PDMP: simulated ', num2str(length(reactions)), ...
                    ' reactions']
            else
                [Tpdmp, Rpdmp, Zpdmp] = simulate_PDMP(z0, t0, t1, ...
                                        prop_func_pdmp, ...
                                        stoch_matrix_pdmp, evo_func_pdmp);
                Xpdmp = repmat((N.^alpha)', size(Zpdmp,1), 1) .* Zpdmp;
                ['PDMP: simulated ', num2str(sum(Rpdmp > 0)), ' reactions']
            end
        end

        figure(2);
        ax(1) = subplot(2, 1, 1);
        if runSSA
            plot_trajectories(Tssa, Xssa, species_str, plot_scaling, ...
                                'SSA Trajectories');
        end
        ax(2) = subplot(2, 1, 2);
        if runPDMP
            plot_trajectories(Tpdmp, Xpdmp, species_str, plot_scaling, ...
                                'PDMP Trajectories');
        end
        if runSSA && runPDMP
            linkaxes([ax(2) ax(1)], 'xy');
        end
    end

    if simulateTrajectoryDistributions
        Mssa = 100;
        Mpdmp = 100;
        dt = 0.1;

        if runSSA
            Tssa = t0:dt:t1;
            Xssa = zeros(Mssa,floor((t1-t0)/dt)+1,s0);
        end
        if runPDMP
            Tpdmp = t0:dt:t1;
            Xpdmp = zeros(Mpdmp,floor((t1-t0)/dt)+1,s0);
        end

        % simulate trajectories
        if runSSA
            'SSA'
            for i = 1:Mssa
                [T, R, X] = simulate_SSA(x0, t0, t1, prop_func, ...
                                         stoch_matrix);
                Y = zeros(floor((t1-t0)/dt)+1,s0);
                for j = 2:length(T)
                    k1 = floor((T(j-1)-t0)/dt)+1;
                    k2 = floor((T(j)-t0)/dt)+1;
                    k2 = min(k2, size(Y,1));
                    k1 = min(k1, k2);
                    Y(k1:k2,:) = repmat(X(j-1,:),k2-k1+1,1);
                end
                Xssa(i,:,:) = Y;
            end
        end
        if runPDMP
            'PDMP'
            a1 = cputime;
            for i = 1:Mpdmp
                if runWithJava
                    [T, X, reactionTimes, reactions] ...
                        = simulate_PDMP_java( ...
                            x0, t0, t1, ...
                            production_stoch_matrix, ...
                            consumption_stoch_matrix, ...
                            prop_rates, ...
                            N, gamma, alpha, beta, ...
                            length(Tpdmp));
                    Xpdmp(i,:,:) = X;
                else
                    [T, R, Z] = simulate_PDMP(z0, t0, t1, ...
                                            prop_func_pdmp, ...
                                            stoch_matrix_pdmp, evo_func_pdmp);
                    X = repmat((N.^alpha)', size(Z,1), 1) .* Z;
                    Y = zeros(floor((t1-t0)/dt)+1,s0);
                    for j = 2:length(T)
                        k1 = floor((T(j-1)-t0)/dt)+1;
                        k2 = floor((T(j)-t0)/dt)+1;
                        k2 = min(k2, size(Y,1));
                        k1 = min(k1, k2);
                        Y(k1:k2,:) = repmat(X(j-1,:),k2-k1+1,1);
                    end
                    Xpdmp(i,:,:) = Y;
                end
            end
            a2 = cputime;
            a2 - a1
        end
        if runSSA
            Xssa_mean = squeeze(mean(Xssa, 1));
            Xssa_std = squeeze(std(Xssa, 1));
        end
        if runPDMP
            Xpdmp_mean = squeeze(mean(Xpdmp, 1));
            Xpdmp_std = squeeze(std(Xpdmp, 1));
        end

        figure(2);
        ax(1) = subplot(2, 1, 1);
        if runSSA
            plot_trajectory_distributions(Tssa, Xssa_mean, Xssa_std, ...
                        species_str, plot_scaling, 'SSA Trajectories');
        end
        ax(2) = subplot(2, 1, 2);
        if runPDMP
            plot_trajectory_distributions(Tpdmp, Xpdmp_mean, Xpdmp_std, ...
                        species_str, plot_scaling, 'PDMP Trajectories');
        end
        if runSSA && runPDMP
            linkaxes([ax(2) ax(1)], 'xy');
        end
    end

end

function plot_trajectories(T, X, species_str, plot_scaling, plot_title)

    % plot trajectory
    cla;
    co = get(gca,'ColorOrder');
    str = {};
    hold on;
    for i = 1:length(X(1,:))
        %plot(T, X(:,i) .* plot_scaling(i), 'Color', co(i, :));
        stairs(T, X(:,i) .* plot_scaling(i), 'Color', co(i, :));
        %str{i} = ['Species ', num2str(i)];
    end
    %legend(str);
    legend(species_str);
    title(plot_title);
    hold off;
    %saveas(gcf(), 'homework3_2A.pdf');

    % Estimation of first and second moment
    dt = T(2:end) - T(1:end-1);
    DT = repmat(dt', [1, size(X, 2)]);
    estimated_first_moment = sum(X(1:end-1,:) .* DT) ./ sum(dt)
    estimated_second_moment = sum(X(1:end-1,:).^2 .* DT) ./ sum(dt)
    estimated_variance = estimated_second_moment ...
                            - estimated_first_moment.^2

end

function plot_trajectory_distributions(T, X_mean, X_std, ...
                        species_str, plot_scaling, plot_title)

    % plot trajectory
    cla;
    co = get(gca,'ColorOrder');
    str = {};
    hold on;
    for i = 1:length(X_mean(1,:))
        %plot(T, X(:,i) .* plot_scaling(i), 'Color', co(i, :));
        plot(T, X_mean(:,i) .* plot_scaling(i), 'Color', co(i, :));
        %str{i} = ['Species ', num2str(i)];
    end
    %legend(str);
    legend(species_str);
    for i = 1:length(X_mean(1,:))
        plot(T, (X_mean(:,i)+X_std(:,i)) .* plot_scaling(i), ':', 'Color', co(i, :));
        plot(T, (X_mean(:,i)-X_std(:,i)) .* plot_scaling(i), ':', 'Color', co(i, :));
    end
    title(plot_title);
    hold off;
    %saveas(gcf(), 'homework3_2A.pdf');

end

function plot_distribution(t1, X, species_str, plot_scaling, plot_title)

    % plot trajectory
    cla;
    hist(X .* repmat(plot_scaling, size(X,1), 1), 1:max(max(X .* repmat(plot_scaling, size(X,1), 1))));
    title(strcat(['Histogram of X(', num2str(t1), ')']));
    legend(species_str);
    title(plot_title);
    hold off;
    %saveas(gcf(), 'homework3_2A.pdf');

    % Estimation of first and second moment
    estimated_first_moment = sum(X) ./ length(X)
    estimated_second_moment = sum(X).^2 ./ length(X)
    estimated_variance = estimated_second_moment ...
                            - estimated_first_moment.^2

end
