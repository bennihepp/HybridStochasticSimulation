function vilar_oscillator_ode()

    initJava();

    import ch.ethz.khammash.hybridstochasticsimulation.examples.*;
    import ch.ethz.khammash.hybridstochasticsimulation.matlab.*;
    import ch.ethz.khammash.hybridstochasticsimulation.models.*;
    import ch.ethz.khammash.hybridstochasticsimulation.networks.*;

    name = 'Vilar Oscillator';
    conf = ExampleConfigurationFactory.getInstance().createExampleConfiguration(name);

    conf.t1 = 100;
    conf.deterministicReactions = [ 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15 ];

    hrn = HybridReactionNetwork(conf.net, conf.deterministicReactions);
    model = HybridReactionNetworkModel(hrn);

    prod = hrn.getProductionStochiometries()';
    cons = hrn.getConsumptionStochiometries()';
    rates = hrn.getRateParameters();

    A = zeros(size(prod,1),size(prod,1));
    B = zeros(size(prod,1),1);

    for r=1:size(prod,2)
        for sp=1:size(A,1)
            if prod(sp,r) > 0
                sc_count = 0;
                for sc=1:size(A,1)
                    if cons(sc,r) > 0
                        sc_count = sc_count + 1;
                    end
                end
                if sc_count >= 2
                    for sc=1:size(A,1)
                        if cons(sc,r) > 0
                            A(sp,sc) = A(sp,sc) + 0.5*rates(r);
                        end
                    end
                elseif sc_count >= 1
                    for sc=1:size(A,1)
                        if cons(sc,r) > 0
                            A(sp,sc) = A(sp,sc) + rates(r);
                        end
                    end
                else
                	B(sp) = B(sp) + rates(r);
                end
            end
        end
    end

    t0 = conf.t0;
    t1 = conf.t1;
    x0 = conf.x0;

    mi = MatlabInterface(model);
    function [dx] = vectorfield(t, x)
        dx = mi.computeVectorField(t, x);
    end

    tic;
    [T, X] = ode15s(@vectorfield, [t0, t1], x0);
    q = toc;
    ['Total execution time: ', num2str(round(q * 1000))]

    plotScales = ones(1,size(X,2));
    Y = X .* repmat(plotScales, size(X,1), 1);
    figure();
    plot(T, Y);

end

function initJava()
    javaHSSBinPath = '../bin/';
    javaHSSLibPath = '../lib/';
    javaaddpath(javaHSSBinPath);
    javaaddpath([javaHSSLibPath, 'commons-lang3-3.1.jar']);
    javaaddpath([javaHSSLibPath, 'commons-math3-3.2.jar']);
    javaaddpath([javaHSSLibPath, 'guava-14.0.1.jar']);
    javaaddpath([javaHSSLibPath, 'JaCoP-3.2.jar']);
    javaaddpath([javaHSSLibPath, 'jamtio.jar']);
    javaaddpath([javaHSSLibPath, 'jcommon-1.0.17.jar']);
    javaaddpath([javaHSSLibPath, 'jfreechart-1.0.14.jar']);
    javaaddpath([javaHSSLibPath, 'jgrapht-jdk1.6.jar']);
    javaaddpath([javaHSSLibPath, 'matlabcontrol-4.1.0.jar']);
end
