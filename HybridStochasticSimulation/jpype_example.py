import os
import platform
import numpy as np
import matplotlib.pyplot as plt
from jpype import *

if __name__ == '__main__':
    classpaths = ["commons-math3-3.2.jar", "jcommon-1.0.17.jar",
                  "JaCoP-3.2.jar", "jfreechart-1.0.14.jar",
                  "jgrapht-jdk1.6.jar", "bin"]
    if platform.system() == "Linux":
        classpath_arg = "-Djava.class.path={}".format(
            ":".join(classpaths))
        startJVM(
            getDefaultJVMPath(),
            "-ea", classpath_arg
        )
    else:
        classpath_arg = "-Djava.class.path={}".format(
            ";".join(classpaths))
        startJVM(
            r"C:\Program Files (x86)\Java\jre7\bin\client\jvm.dll",
            "-ea", classpath_arg
        )
    org = JPackage("org")
    int = org.apache.commons.math3.ode.nonstiff.DormandPrince853Integrator(JDouble(1.0),JDouble(1.0),JDouble(1.0),JDouble(1.0))

    ReactionNetwork = JClass("ReactionNetwork")
    HybridReactionNetwork = JClass("HybridReactionNetwork")
    HybridReactionNetworkModel = JClass("HybridReactionNetworkModel")
    PDMPModel = JClass("PDMPModel")
    PDMPModelSimulator = JClass("PDMPModelSimulator")
    PDMPOutputModel = JClass("PDMPOutputModel")

    net = ReactionNetwork(4, 2)
    net.setStochiometry(0, 0, 0, 2);
    net.setStochiometry(1, 0, 1, 0);
    net.setStochiometry(2, 0, 0, 0);
    net.setStochiometry(3, 0, 0, 0);
    net.setStochiometry(0, 1, 0, 1);
    net.setStochiometry(1, 1, 0, 0);
    net.setStochiometry(2, 1, 0, 1);
    net.setStochiometry(3, 1, 1, 0);
    net.setRateParameter(0, 1e-7);
    net.setRateParameter(1, 1e-7);
    N = 1e6
    gamma = 0
    alpha = [1, 1, 0, 0]
    beta = [-1, -1]
    hrn = HybridReactionNetwork(net, N, gamma, alpha, beta)
    hrnModel = HybridReactionNetworkModel(hrn)
    pdmpModel = PDMPModel(hrnModel, hrnModel)
    pdmpSolver = PDMPModelSimulator()
    cm = PDMPOutputModel()

    pdmpSolver.addStepHandler(cm)
    pdmpSolver.addReactionHandler(cm)
    t0 = 0.0
    t1 = 100.0
    x0 = [1.0e6, 0.0, 10.0, 0.0]
    z0 = hrn.scaleState(x0)
    z1 = list(z0)
    K = 100
    M = 1000
    T = np.linspace(t0, t1, M)
    X = np.zeros((K, T.shape[0], len(z0)))
    for k in xrange(K):
        pdmpSolver.solve(pdmpModel, t0, z0, t1, z1)
        for i in xrange(T.shape[0]):
            cm.setInterpolatedTime(T[i])
            z = cm.getInterpolatedState()
            x = hrn.recoverState(z)
            X[k,i,:] += x
    pdmpSolver.removeStepHandler(cm)
    pdmpSolver.removeReactionHandler(cm)

    Xmean = X.mean(0)
    XstdDev = X.std(0)

    plotScale = np.array([1.0e-5, 1.0e-5, 1.0, 1.0])
    ax = plt.axes()
    ax.set_color_cycle(["b", "g", "r", "c"])
    plt.plot(T, Xmean * plotScale)
    plt.plot(T, (Xmean + XstdDev) * plotScale, ":")
    plt.plot(T, (Xmean - XstdDev) * plotScale, ":")
    plt.legend(["S1", "S2", "S3", "S4"])
    plt.show()

    shutdownJVM()
