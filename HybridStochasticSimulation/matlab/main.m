function main()

    initJava();

    import ch.ethz.khammash.hybridstochasticsimulation.examples.*;
    import ch.ethz.khammash.hybridstochasticsimulation.matlab.*;
    import ch.ethz.khammash.hybridstochasticsimulation.models.*;
    import ch.ethz.khammash.hybridstochasticsimulation.networks.*;

    name = 'Heat Shock Response';
    conf = ExampleConfigurationFactory.getInstance().createExampleConfiguration(name);

%     conf.N = 10000;
%     conf.epsilon = 0.1;
%     conf.xi = 0.5;
%     conf.delta = 0.5;
%     conf.gamma = 0;
    conf.N = 100;
    conf.epsilon = 0.5;
    conf.xi = 0.1;
    conf.delta = 0.1;
    conf.gamma = 2;
    conf.t1 = 100000;
    conf.t1 = 10000000;
    conf.t1 = 2*10000;

    hrn = AdaptiveMSHRN(conf.net, conf.N, conf.gamma, conf.alpha, conf.beta);
    hrn.setDelta(conf.delta);
    hrn.setEpsilon(conf.epsilon);
    hrn.setXi(conf.xi);
    hrn.setTolerance(conf.tolerance);
    model = AdaptiveMSHRNModel(hrn);

    numOfTimePoints = 1001;

    t0 = conf.t0;
    t1 = conf.t1;
    x0 = conf.x0;
    tic;
    [T, X] = PDMPSimulator(model, t0, x0, t1, numOfTimePoints);
    q = toc;
    ['Total execution time: ', num2str(round(q * 1000))]

    plotScales = ones(size(X,1),1);
    Y = X .* repmat(plotScales, 1, numOfTimePoints);
    figure();
    plot(T, Y);

    figure();
%     subplot(3, 1, 1);
%     plot(T, X(2,:));
%     subplot(3, 1, 2);
%     plot(T, X(3,:));
%     subplot(3, 1, 3);
%     plot(T, X(8,:));
    rows = ceil(sqrt(size(X,1)));
    for s = 1:size(X,1)
        subplot(rows, rows, s);
        plot(T, X(s,:));
        title(['S', int2str(s)]);
    end

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
