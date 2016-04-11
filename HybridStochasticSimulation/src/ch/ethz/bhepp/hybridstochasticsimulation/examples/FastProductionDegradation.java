package ch.ethz.bhepp.hybridstochasticsimulation.examples;

import java.util.Arrays;

import org.apache.commons.math3.random.MersenneTwister;
import org.apache.commons.math3.random.RandomDataGenerator;

import ch.ethz.bhepp.hybridstochasticsimulation.networks.DefaultUnaryBinaryReactionNetwork;


// Fast expression and degradation example.
// species:
//  S1
//  S2
//  P1
//  P2
// reactions:
//  R0:  S1 -> S2       [0.002]
//  R1:  S2 -> S1       [0.002]
//  R2:  S1 -> S1 + P1  [20.0]
//  R3:  S2 -> S2 + P2  [20.0]
//  R4:  P1 -> ~        [0.002]
//  R5:  P2 -> ~        [0.002]

public class FastProductionDegradation extends SimulationConfiguration {

    public FastProductionDegradation() {

        int[][] productionStochiometries = {
        //    S1  S2  P1  P2
            { 0,  1,  0,  0 }, // R0
            { 1,  0,  0,  0 }, // R1
            { 1,  0,  1,  0 }, // R2
            { 0,  1,  0,  1 }, // R3
            { 0,  0,  0,  0 }, // R4
            { 0,  0,  0,  0 }, // R5
        };

        int[][] consumptionStochiometries = {
        //    S1  S2  P1  P2
            { 1,  0,  0,  0 }, // R0
            { 0,  1,  0,  0 }, // R1
            { 1,  0,  0,  0 }, // R2
            { 0,  1,  0,  0 }, // R3
            { 0,  0,  1,  0 }, // R4
            { 0,  0,  0,  1 }, // R5
        };

        double[] x0 = { 500, 0, 0, 0 };

        double[] rateParameters = new double[productionStochiometries.length];
        /* Reaction propensities */
        rateParameters[0] = 0.002;
        rateParameters[1] = 0.002;
        rateParameters[2] = 20.0;
        rateParameters[3] = 20.0;
        rateParameters[4] = 0.002;
        rateParameters[5] = 0.002;

        double t0 = 0.0;
        double t1 = 1000.0;
        double N = 100;
        double gamma = 0;
        double[] alpha = new double[x0.length];
        double[] beta = new double[rateParameters.length];

        DefaultUnaryBinaryReactionNetwork net = new DefaultUnaryBinaryReactionNetwork(x0.length, rateParameters.length);
        net.setStoichiometries(productionStochiometries, consumptionStochiometries);
        net.setRateParameters(rateParameters);
        String[] speciesNames = {
            "S1",
            "S2",
            "P1",
            "P2",
        };
        double[] plotScales = new double[x0.length];
        Arrays.fill(plotScales, 1.0);
        int[] importantSpecies = { 0, 1 };

        this.net = net;
        this.N = N;
        this.gamma = gamma;
        this.alpha = alpha;
        this.beta = beta;
        this.t0 = t0;
        this.tf = t1;
        this.x0 = x0;
        this.plotScales = plotScales;
        this.speciesNames = speciesNames;
        this.importantSpecies = importantSpecies;
        this.rng = new MersenneTwister();
        this.rdg = new RandomDataGenerator(this.rng);
    }

}
