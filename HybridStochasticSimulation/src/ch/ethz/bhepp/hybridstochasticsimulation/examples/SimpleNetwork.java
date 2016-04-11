package ch.ethz.bhepp.hybridstochasticsimulation.examples;

import java.util.Arrays;

import org.apache.commons.math3.random.MersenneTwister;
import org.apache.commons.math3.random.RandomDataGenerator;

import ch.ethz.bhepp.hybridstochasticsimulation.networks.DefaultUnaryBinaryReactionNetwork;

//Species:
//S0: A
//S1: B
//Reactions:
//R0: ~ -> A
//R1: A -> ~
//R2: A -> A + B
//R3: B -> ~

public class SimpleNetwork extends SimulationConfiguration {

	public SimpleNetwork() {

        int[][] productionStochiometries = {
        //    S0  S1
            { 1,  0 }, // R0
            { 0,  0 }, // R1
            { 1,  1 }, // R2
            { 0,  0 }, // R3
        };

        int[][] consumptionStochiometries = {
        //    S0  S1
            { 0,  0 }, // R0
            { 1,  0 }, // R1
            { 1,  0 }, // R2
            { 0,  1 }, // R3
        };

        double[] x0 = { 0, 0 };

        double[] rateParameters = new double[productionStochiometries.length];
        /* Reaction propensities */
        rateParameters[0] = 0.1;
        rateParameters[1] = 0.005;
        rateParameters[2] = 10.0;
        rateParameters[3] = 0.1;

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
            "A",
            "B",
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
