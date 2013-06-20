function [T, X] = PDMPSimulator( ...
                        model, t0, x0, t1, numOfTimePoints)

    import ch.ethz.khammash.hybridstochasticsimulation.matlab.*;

    model.initialize(t0, x0);
    mi = MatlabInterface(model);

    % set initial time and state
    t = t0;
    z0 = model.getNetwork().scaleState(x0);
    x = z0;
    x = vertcat(x, [0; 0]);

    % T will hold time points, X will hold state values
    % and R will hold reaction indices
    T = linspace(t0, t1, numOfTimePoints);
%     R = zeros(size(T));
    X = zeros(length(x0),length(T));

    nextOutputIndex = 1;
    function [status] = ode_output_function(t, x, opt)
        switch opt
            case 'init'
            case 'done'
            otherwise
                for i = 1:length(t)
                    if t(i) ~= T(nextOutputIndex)
                        break;
                    end
                    z = x(1:end-2,i);
                    recoveredX = mi.recoverStateVector(t(i), z);
                    T(nextOutputIndex) = t(i);
                    X(:,nextOutputIndex) = recoveredX;
                    nextOutputIndex = nextOutputIndex + 1;
                end
%                 if t(1) == T(nextOutputIndex)
%                     T(nextOutputIndex:nextOutputIndex+length(t)-1) = t;
%                     X(:,nextOutputIndex:nextOutputIndex+length(t)-1) = x(1:end-2);
%                     nextOutputIndex = nextOutputIndex + 1;
%                 end
                status = 0;
        end
    end
    ode_output_function(t, x, '');

    function [dx] = vectorfield(t, x)
        if t > 11691.0875
            t = t;
        end
        dx = mi.computeVectorField(t, x);
    end

%     lowerBounds = -1 * ones(length(z0), 1);
%     upperBounds = 100 * ones(length(z0), 1);
%     bounds = [lowerBounds, upperBounds];
    numOfOptionalObservers = mi.getNumberOfOptionalObservers();
    function [values, isterminal, direction] = ode_events(t, x)
        isterminal = ones(1, 1 + numOfOptionalObservers);
        direction = zeros(1, 1 + numOfOptionalObservers);
        value1 = x(end) - x(end - 1);
        value2 = mi.computeOptionalObserverValues(t, x);
        values = [value1, value2'];
%         values = zeros(length(z0) + 1);
%         for s=1:length(z0)
%             lowerBound = bounds(s, 1);
%             upperBound = bounds(s, 2);
%             if lowerBound < 0
%                 values(s) = x(s) - upperBound;
%             else
%                 values(s) = (x(s) - lowerBound) * (x(s) - upperBound);
%             end
%         end
%         values(end) = x(end) - x(end - 1);
    end

    ode_options = odeset('RelTol', 1e-6, 'AbsTol', 1e-6, 'Events', @ode_events, 'OutputFcn', @ode_output_function);

%     T(1) = t0;
%     R(1) = 0;
%     X(1,:) = x0;

%     % index variable i
%     i = 1;
    % iteration variable
    j = 1;

    tic;
    while t < t1
        x(end - 1) = 0;
        x(end) = -log(rand(1));
        % evolve deterministic part until first discrete reaction
        %tspan = [t, t1];
        while t < t1
            x = mi.handleOptionalEvent(t, x);
            tspan = T(nextOutputIndex:end);
            if tspan(1) > t
                tspan = [t, tspan];
            end
    %         tspan = T(i:end);
    %         if T(i) ~= t
    %             tspan = [t, tspan];
    %         end
            [TT, XX, TE, XE, IE] = ode15s(@vectorfield, tspan, x, ode_options);
    %         if T(i) ~= t
    %             TT = TT(2:end);
    %             XX = XX(2:end,:);
    %         end
    %         T(i:end) = TT(:);
    %         X(i:end,:) = XX(:,:);
    %         i = i + length(TT);
            t = TT(end);
            x = XX(end,:);
            if length(IE) > 0
                t = t;
            end
            if length(IE) > 0
                if IE(1) == 1
                    break;
                else
                    mi.reportOptionalEvent(0, t, x);
                end
            end
        end

        if t >= t1
            break;
        end

        p = mi.computeTransitionMeasure(t, x);
        psum = sum(p);
        u = rand(1);
        w = 0;
        r = 0;
        for l = 1:length(p)
            w = w + p(l) / psum;
            if u < w
                r = l;
                x = mi.updateState(r - 1, t, x);
                break;
            end
        end
        if r > 0
            x = mi.checkAndHandleOptionalEvent(t, x);
        end
        j = j + 1;
        if mod(j, 100)
%             [j, t, 100 * t / t1]
%             j
%             t
%             dt = 100 * (t-t0) / (t1-t0);
%             dt
        end
    end
    q = toc;
    ['Execution time: ', num2str(round(q * 1000))]

end
