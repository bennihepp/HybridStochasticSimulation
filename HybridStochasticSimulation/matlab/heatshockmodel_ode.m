function heatshockmodel_ode()

    N0 = 100;
    gamma = 1;
    timescale = N0^gamma;

    [s0, r0, stoch_matrix, prop_func, x0] = heatshockmodel();
    t0 = 0;
    t1 = 10000;
    tspan = [t0; t1] / timescale;
    z0 = x0 / timescale;

    function xdot = ode(t, x)
        xdot = zeros(size(x));
        props = prop_func(x * timescale);
        for r=1:r0
            prop = props(r);
            for s=1:s0
                xdot(s) = xdot(s) + prop * stoch_matrix(r, s);
            end
            xdot(s) = xdot(s) / timescale;
        end
    end

    options = odeset('RelTol',1e-3,'AbsTol',1e-6,'MaxStep',0.1*abs(t1 - t0));
    tic;
    [TAU, Z] = ode23s(@ode, tspan, z0, options);
    X = timescale * Z;
    T = timescale * TAU;
    q = toc;
    ['Total execution time: ', num2str(round(q * 1000))]

    figure();
    plot(T, X);
    title(['\gamma = ', num2str(gamma)]);

end

function [s0, r0, stoch_matrix, prop_func, x0] = heatshockmodel()
    % model parameters

    production_stoch_matrix = [ ...
        0,  0,  0,  0,  0,  0,  0,  1,  0; ... % R1
        0,  0,  1,  0,  0,  0,  0,  0,  0; ... % R2
        0,  1,  0,  0,  0,  0,  0,  0,  0; ... % R3
        1,  1,  0,  0,  0,  0,  0,  0,  0; ... % R4
        0,  1,  0,  0,  1,  0,  0,  0,  0; ... % R5
        0,  1,  0,  1,  0,  0,  0,  0,  0; ... % R6
        0,  1,  0,  0,  0,  1,  0,  0,  0; ... % R7
        0,  1,  0,  0,  0,  1,  0,  0,  0; ... % R8
        0,  0,  0,  0,  0,  0,  1,  0,  0; ... % R9
        0,  0,  0,  0,  0,  0,  0,  0,  1; ... % R10
        0,  0,  0,  0,  0,  0,  0,  0,  0; ... % R11
        0,  0,  0,  0,  0,  1,  0,  1,  0; ... % R12
        1,  0,  0,  0,  0,  0,  0,  0,  0; ... % R13
        0,  0,  0,  0,  0,  0,  0,  0,  0; ... % R14
        0,  0,  0,  1,  0,  1,  0,  0,  0; ... % R15
        0,  0,  0,  0,  0,  0,  0,  0,  0; ... % R16
        0,  0,  0,  0,  0,  0,  0,  0,  0; ... % R17
        0,  0,  0,  0,  0,  0,  0,  0,  0; ... % R18
    ];
    consumption_stoch_matrix = [ ...
        0,  0,  0,  0,  0,  0,  0,  0,  0; ... % R1
        0,  1,  0,  0,  0,  0,  0,  0,  0; ... % R2
        0,  0,  1,  0,  0,  0,  0,  0,  0; ... % R3
        1,  0,  0,  0,  0,  0,  0,  0,  0; ... % R4
        0,  0,  1,  0,  0,  0,  0,  0,  0; ... % R5
        0,  0,  1,  0,  0,  0,  0,  0,  0; ... % R6
        0,  0,  1,  0,  0,  0,  0,  0,  0; ... % R7
        0,  0,  0,  0,  0,  0,  1,  0,  0; ... % R8
        0,  1,  0,  0,  0,  1,  0,  0,  0; ... % R9
        0,  0,  0,  0,  0,  1,  0,  1,  0; ... % R10
        0,  0,  0,  0,  0,  0,  0,  1,  0; ... % R11
        0,  0,  0,  0,  0,  0,  0,  0,  1; ... % R12
        0,  0,  0,  0,  0,  0,  0,  0,  0; ... % R13
        1,  0,  0,  0,  0,  0,  0,  0,  0; ... % R14
        0,  0,  0,  1,  0,  0,  1,  0,  0; ... % R15
        0,  0,  0,  0,  1,  0,  0,  0,  0; ... % R16
        0,  0,  0,  0,  0,  1,  0,  0,  0; ... % R17
        0,  0,  0,  1,  0,  0,  0,  0,  0; ... % R18
    ];
    stoch_matrix = production_stoch_matrix - consumption_stoch_matrix;
    r0 = size(stoch_matrix,1);
    s0 = size(stoch_matrix,2);

%     x0 = zeros(s0, 1);
    x0 = [ 10; 1; 1; 93; 172; 54; 7; 50; 0 ];

    rates = [ 4.0, 0.7, 0.13, 0.007, 0.0063, 0.00488, 0.00488, 4.4E-4, 3.62E-4, 3.62E-4, 9.99E-5, 4.4E-5, 1.4E-5, 1.4E-6, 1.42E-6, 1.8E-8, 6.4E-10, 7.4E-11 ]';
    choiceIndices = [
        -1, -1;
         2, -1;
         3, -1;
         1, -1;
         3, -1;
         3, -1;
         3, -1;
         7, -1;
         2,  6;
         6,  8;
         8, -1;
         9, -1;
        -1, -1;
         1, -1;
         7, -1;
         5, -1;
         6, -1;
         4, -1;
    ];
    function props = prop_function(x)
        props = zeros(size(rates));
        for r=1:r0
            choiceIndex1 = choiceIndices(r, 1);
            choiceIndex2 = choiceIndices(r, 2);
            v = rates(r);
            if choiceIndex1 > 0
                v = v * x(choiceIndex1);
            end
            if choiceIndex2 > 0
                v = v * x(choiceIndex2);
            end
            props(r) = v;
        end
    end
    prop_func = @prop_function;
%     prop_func = @(x) [
%         rates(1); ...
%         rates(2)*x(2); ...
%         rates(3)*x(3); ...
%         rates(4)*x(1); ...
%         rates(5)*x(3); ...
%         rates(6)*x(3); ...
%         rates(7)*x(3); ...
%         rates(8)*x(7); ...
%         rates(9)*x(2)*x(6); ...
%         rates(10)*x(6)*x(8); ...
%         rates(11)*x(8); ...
%         rates(12)*x(9); ...
%         rates(13); ...
%         rates(14)*x(1); ...
%         rates(15)*x(7); ...
%         rates(16)*x(5); ...
%         rates(17)*x(6); ...
%         rates(18)*x(4);
%     ];

end
