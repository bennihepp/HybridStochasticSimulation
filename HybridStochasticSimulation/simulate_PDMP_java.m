function [T, X, reactionTimes, reactions] = simulate_PDMP_java( ...
                        x0, t0, t1, ...
                        production_stoch_matrix, ...
                        consumption_stoch_matrix, ...
                        prop_rates, ...
                        N, gamma, alpha, beta, ...
                        steps)

    r0 = size(production_stoch_matrix,1);
    s0 = size(production_stoch_matrix,2);

    mi = MatlabInterface();

    rn = mi.createReactionNetwork(s0, r0);
    rn.setStochiometries( ...
        production_stoch_matrix, consumption_stoch_matrix);
    rn.setRateParameters(prop_rates);

    hrn = mi.createHybridReactionNetwork(rn, N, gamma, alpha, beta);

    hrnModel = mi.createHybridReactionNetworkModel(hrn);
    pdmpModel = mi.createPDMPModel(hrnModel, hrnModel);
    pdmpModelSimulator = mi.createPDMPModelSimulator();
    cm = mi.createPDMPContinuousOutputModel();

    z0 = hrn.scaleState(x0);
    z1 = zeros(size(z0));

    pdmpModelSimulator.addStepHandler(cm);
    pdmpModelSimulator.addReactionHandler(cm);
    pdmpModelSimulator.simulate(pdmpModel, t0, z0, t1, z1);
    pdmpModelSimulator.removeReactionHandler(cm);
    pdmpModelSimulator.removeStepHandler(cm);

    T = linspace(t0, t1, steps);
    X = zeros(length(T), length(z0));
    for i = 1:length(T)
        cm.setInterpolatedTime(T(i));
        z = cm.getInterpolatedState();
        x = hrn.recoverState(z);
        X(i,:) = x;
    end

    reactionTimes = zeros(cm.getNumberOfReactions(), 1);
    reactions = zeros(cm.getNumberOfReactions(), 1);
    for i = 1:cm.getNumberOfReactions()
        reactions(i) = cm.getIndexOfReaction(i-1);
        reactionTimes(i) = cm.getTimeOfReaction(i-1);
    end

end

function [value, isterminal, direction] = ode_events(t, x)
    value = x(end) - x(end - 1);
    isterminal = 1;
    direction = 0;
end
