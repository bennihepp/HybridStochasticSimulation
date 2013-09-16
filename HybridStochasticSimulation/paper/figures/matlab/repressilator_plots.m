addpath('../../../matlab/hdf5');


indices = [ 1, 3, 5 ];


[Tstoch, Xstoch] = readSimulationData( ...
    '../../../jobs/outputs/Repressilator_stoch_10000_full.h5', ...
    'simulations' ...
);

figure;
plot_single_trace(Tstoch, Xstoch, 1, indices);
title('Repressilator stochastic');
legend('m1', 'm2', 'm3');
save_plot('../repressilator_stochastic_single');
figure;
plot_mean(Tstoch, Xstoch, indices);
title('Repressilator stochastic');
legend('m1', 'm2', 'm3');
save_plot('../repressilator_stochastic_mean');


[Tadaptive, Xadaptive] = readSimulationData( ...
    '../../../jobs/outputs/Repressilator_adaptive_10000_full.h5', ...
    'simulations' ...
);

figure;
plot_single_trace(Tadaptive, Xadaptive, 1, indices);
title('Repressilator adaptive');
legend('m1', 'm2', 'm3');
save_plot('../repressilator_adaptive_single');
figure;
plot_mean(Tadaptive, Xadaptive, indices);
title('Repressilator adaptive');
legend('m1', 'm2', 'm3');
save_plot('../repressilator_adaptive_mean');


[Tmspdmp, Xmspdmp] = readSimulationData( ...
    '../../../jobs/outputs/Repressilator_mspdmp_10000_full.h5', ...
    'simulations' ...
);

figure;
plot_single_trace(Tmspdmp, Xmspdmp, 1, indices);
title('Repressilator MSPDMP');
legend('m1', 'm2', 'm3');
save_plot('../repressilator_mspdmp_single');
figure;
plot_mean(Tmspdmp, Xmspdmp, indices);
title('Repressilator MSPDMP');
legend('m1', 'm2', 'm3');
save_plot('../repressilator_mspdmp_mean');
