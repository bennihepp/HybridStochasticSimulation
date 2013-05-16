% model parameters
k1 = 10;
k2 = 100;
k3 = 50;
k4 = 5;
x0 = [0; 0];

stoch_matrix = [1, 0; ...
                -1, 1; ...
                1, -1; ...
                0, -1];
stoch_size = size(stoch_matrix);
r0 = stoch_size(1);
s0 = stoch_size(2);

prop_rates = [k1; k2; k3; k4];
prop_func = @(x) [prop_rates(1); prop_rates(2)*x(1); ...
                  prop_rates(3)*x(2); prop_rates(4)*x(2)];
evo_func = @(t, x) [-prop_rates(1) - prop_rates(2)*x(1) + prop_rates(3)*x(2); ...
                    prop_rates(2)*x(1) - prop_rates(3)*x(2) - prop_rates(4)*x(2)];
species_str = {'A', 'B'};
plot_scaling = [1, 1];
