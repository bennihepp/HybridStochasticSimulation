package ch.ethz.bhepp.hybridstochasticsimulation.examples;

import java.util.Arrays;

import org.apache.commons.math3.random.MersenneTwister;
import org.apache.commons.math3.random.RandomDataGenerator;

import ch.ethz.bhepp.hybridstochasticsimulation.networks.DefaultUnaryBinaryReactionNetwork;


// Fast expression and degradation example.
// species:
//  S1
// reactions:
//  R0:  ~  -> S1  [100]
//  R1:  S1 -> ~   [1]

public class VeryFastProductionDegradation extends SimulationConfiguration {

    public VeryFastProductionDegradation() {

        int[][] productionStochiometries = {
        //    S1
            { 1 }, // R0
            { 0 }, // R1
        };

        int[][] consumptionStochiometries = {
        //    S1
            { 0 }, // R0
            { 1 }, // R1
        };

        double[] x0 = { 0 };

        double[] rateParameters = new double[productionStochiometries.length];
        /* Reaction propensities */
        rateParameters[0] = 100;
        rateParameters[1] = 1;

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
