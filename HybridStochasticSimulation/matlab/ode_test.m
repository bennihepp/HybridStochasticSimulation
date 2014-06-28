function ode_test()

    N0 = 100;
    gamma = 1;
    timescale = N0^gamma;

    x0 = [ 1; 0; 0 ];
    t0 = 0;
    t1 = 4*10^5;
    tspan = [t0; t1] / timescale;
%     z0 = x0 / timescale;
    z0 = x0;

    function xdot = ode(t, x)
        xdot = zeros(size(x));
        xdot(1) = -0.03 * x(1) + 0.5 * x(3);
        xdot(2) = -0.3 * x(1) - 0.2 * x(2);
        xdot(3) = 0.02 * x(2) - 0.05 * x(3);
        xdot = xdot * timescale;
%         xdot(1) = -0.04*x(1) + 10^4*x(2)*x(3);
%         xdot(2) = 0.04*x(1) - 10^4*x(2)*x(3) - 3*10^7*x(2)^2;
%         xdot(3) = 3*10^7*x(2)^2;
    end

    options = odeset('RelTol',1e-3,'AbsTol',1e-6,'MaxStep',0.1*abs(t1 - t0));
    tic;
    [TAU, Z] = ode23s(@ode, tspan, z0, options);
    q = toc;
    ['Total execution time: ', num2str(round(q * 1000))]
%     X = timescale * Z;
    X = Z;
    T = timescale * TAU;

    figure();
    plot(T, X);
    title(['\gamma = ', num2str(gamma)]);

end
