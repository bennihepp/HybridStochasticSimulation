% model parameters
kgon = 0.2;
kgoff = 0.1;
kr = 5;
yr = 5;
kp = 2000;
yp = 1;

stoch_matrix = [1, -1, 0, 0; ...
                -1, 1, 0, 0; ...
                0, 0, 1, 0; ...
                0, 0, -1, 0; ...
                0, 0, 0, 1; ...
                0, 0, 0, -1];
% nu_prime_matrix = zeros(size(stoch_matrix));
% nu_prime_matrix(stoch_matrix > 0) = stoch_matrix(stoch_matrix > 0);
% nu_matrix = zeros(size(stoch_matrix));
% nu_matrix(stoch_matrix < 0) = -stoch_matrix(stoch_matrix < 0);

prop_func = @(x) [kgoff*x(2); kgon*x(1); kr*x(2); ...
                    yr*x(3); kp*x(3); yp*x(4)];
evo_func = @(t, x) [kgoff*x(2); kgon*x(1); kr*x(2) - yr*x(3); ...
                    kp*x(3) - yp*x(4)];
species_str = {'GeneOff', 'GeneOn', 'mRNA', 'Protein'};
