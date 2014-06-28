% model parameters

N = 10^6;

k1 = 10^-7;
k2 = 10^-7;
x0 = [10^6; 0; 10; 0];

production_stoch_matrix = [0, 1, 0, 0; ...
                        0, 0, 0, 1];
consumption_stoch_matrix = [2, 0, 0, 0; ...
                        1, 0, 1, 0];
stoch_matrix = production_stoch_matrix - consumption_stoch_matrix;
stoch_size = size(stoch_matrix);
r0 = stoch_size(1);
s0 = stoch_size(2);
% nu_prime_matrix = zeros(size(stoch_matrix));
% nu_prime_matrix(stoch_matrix > 0) = stoch_matrix(stoch_matrix > 0);
% nu_matrix = zeros(size(stoch_matrix));
% nu_matrix(stoch_matrix < 0) = -stoch_matrix(stoch_matrix < 0);

prop_rates = [k1; k2];
prop_func = @(x) [prop_rates(1)*x(1)*(x(1)-1); prop_rates(2)*x(1)*x(3)];
evo_func = @(t, x) [-k1*x(1)^2 - k2*x(1)*x(3); ...
                    k1*x(1)^2; ...
                    -k2*x(1)*x(3); ...
                    k2*x(1)*x(3)];
species_str = {'A', 'B', 'C', 'D'};
plot_scaling = [10^-5, 10^-5, 1, 1];
