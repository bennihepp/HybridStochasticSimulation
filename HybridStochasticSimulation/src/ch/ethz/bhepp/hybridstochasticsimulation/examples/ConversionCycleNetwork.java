package ch.ethz.bhepp.hybridstochasticsimulation.examples;


import org.apache.commons.math3.random.MersenneTwister;
import org.apache.commons.math3.random.RandomDataGenerator;

import ch.ethz.bhepp.hybridstochasticsimulation.networks.DefaultUnaryBinaryReactionNetwork;


public class ConversionCycleNetwork extends SimulationConfiguration {

	public ConversionCycleNetwork() {
		int[] deterministicReactions = { 0, 1 };
		double N = 1e5;
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
		net.setStoichiometries(productionStochiometries, consumptionStochiometries);
		net.setRateParameters(rateParameters);
		String[] speciesNames = { "S1", "S2", "S3" };

		this.net = net;
		this.deterministicReactions = deterministicReactions;
		this.N = N;
		this.gamma = gamma;
		this.alpha = alpha;
		this.beta = beta;
		this.t0 = t0;
		this.tf = t1;
		this.x0 = x0;
		this.plotScales = plotScales;
		this.speciesNames = speciesNames;
		this.rng = new MersenneTwister();
		this.rdg = new RandomDataGenerator(this.rng);
	}

}
