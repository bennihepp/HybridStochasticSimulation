addpath('../../../matlab/hdf5');


indices = [ 4, 5 ];


% [Tstoch, Xstoch] = readSimulationData( ...
%     '../../../jobs/outputs/ToggleSwitch_stoch_10000.h5', ...
%     'simulations' ...
% );
% 
% figure;
% plot_single_trace(Tstoch, Xstoch, 1, indices);
% title('ToggleSwitch stochastic');
% legend('m1', 'm2', 'm3');
% save_plot('../toggleswitch_stochastic_single');
% figure;
% plot_mean(Tstoch, Xstoch, indices);
% title('ToggleSwitch stochastic');
% legend('m1', 'm2', 'm3');
% save_plot('../toggleswitch_stochastic_mean');


[Tadaptive, Xadaptive] = readSimulationData( ...
    '../../../jobs/outputs/ToggleSwitch_adaptive_10000.h5', ...
    'simulations' ...
);

figure;
t = 1000;
plot_timepoint_distribution(t, Tadaptive, Xadaptive, indices, {'m1', 'm2', 'm3'});
title(['ToggleSwitch Adaptive t=', t]);
save_plot('../toggleswitch_adaptive_distribution');
figure;
plot_single_trace(Tadaptive, Xadaptive, 1, indices);
title('ToggleSwitch adaptive');
legend('m1', 'm2', 'm3');
save_plot('../toggleswitch_adaptive_single');
figure;
plot_mean(Tadaptive, Xadaptive, indices);
title('ToggleSwitch adaptive');
legend('m1', 'm2', 'm3');
save_plot('../toggleswitch_adaptive_mean');


[Tmspdmp, Xmspdmp] = readSimulationData( ...
    '../../../jobs/outputs/ToggleSwitch_mspdmp_10000.h5', ...
    'simulations' ...
);

figure;
t = 1000;
t = plot_timepoint_distribution(t, Tmspdmp, Xmspdmp, indices, {'m1', 'm2', 'm3'});
title(['ToggleSwitch MSPDMP t=', t]);
legend('m1', 'm2', 'm3');
save_plot('../toggleswitch_mspdmp_distribution');
figure;
plot_single_trace(Tmspdmp, Xmspdmp, 1, indices);
title('ToggleSwitch MSPDMP');
legend('m1', 'm2', 'm3');
save_plot('../toggleswitch_mspdmp_single');
figure;
plot_mean(Tmspdmp, Xmspdmp, indices);
title('ToggleSwitch MSPDMP');
legend('m1', 'm2', 'm3');
save_plot('../toggleswitch_mspdmp_mean');
