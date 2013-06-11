package ch.ethz.khammash.hybridstochasticsimulation.examples;


import org.apache.commons.math3.random.MersenneTwister;
import org.apache.commons.math3.random.RandomDataGenerator;

import ch.ethz.khammash.hybridstochasticsimulation.networks.DefaultUnaryBinaryReactionNetwork;


public class ConversionCycleNetwork extends ExampleNetwork {

	public ConversionCycleNetwork() {
		int[] continuousSpecies = { 0, 1 };
		double N = 1e5;
		double deltaR = 0.5;
		double deltaS = 0.5;
		double epsilon = 0.5;
		double gamma = 0;
		double[] alpha = { 0, 1, 1 };
		double[] beta = { 0, 0, 0, 0 };
		double t0 = 0.0;
		double t1 = 100.0;
		double[] x0 = { 1, 10000, 10000 };
		double[] plotScales = { 10, 1, 1 };
		int[][] productionStochiometries = {
				{ 0, 1, 0 },
				{ 0, 0, 1 },
				{ 1, 0, 0 },
				{ 0, 1, 0 },
		};
		int[][] consumptionStochiometries = {
				{ 1, 0, 0 },
				{ 0, 1, 0 },
				{ 0, 0, 1 },
				{ 0, 0, 1 },
		};
		double[] rateParameters = {
				10,
				100,
				0.1,
				100,
		};
		DefaultUnaryBinaryReactionNetwork net = new DefaultUnaryBinaryReactionNetwork(productionStochiometries[0].length, productionStochiometries.length);
		net.setStochiometries(productionStochiometries, consumptionStochiometries);
		net.setRateParameters(rateParameters);
		String[] speciesNames = { "S1", "S2", "S3" };

		this.net = net;
		this.continuousSpecies = continuousSpecies;
		this.N = N;
		this.deltaR = deltaR;
		this.deltaS = deltaS;
		this.epsilon = epsilon;
		this.gamma = gamma;
		this.alpha = alpha;
		this.beta = beta;
		this.t0 = t0;
		this.t1 = t1;
		this.x0 = x0;
		this.plotScales = plotScales;
		this.speciesNames = speciesNames;
		this.rng = new MersenneTwister();
		this.rdg = new RandomDataGenerator(this.rng);
	}

}
