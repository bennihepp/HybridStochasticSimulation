package ch.ethz.khammash.hybridstochasticsimulation.examples;

import org.apache.commons.math3.random.MersenneTwister;
import org.apache.commons.math3.random.RandomDataGenerator;

import ch.ethz.khammash.hybridstochasticsimulation.networks.DefaultUnaryBinaryReactionNetwork;

public class TrivialNetwork extends SimulationConfiguration {

	public TrivialNetwork() {
		int[] deterministicReactions = { };
		double N = 1e6;
		double gamma = 0;
		double[] alpha = { 0 };
		double[] beta = { 0, 0 };
		double t0 = 0.0;
		double t1 = 100.0;
		double[] x0 = { 0 };
		double[] plotScales = { 1 };
		int[][] productionStochiometries = {
				{ 1 },
				{ 0 },
		};
		int[][] consumptionStochiometries = {
				{ 0 },
				{ 1 },
		};
		double[] rateParameters = {
				5.0,
				2.0,
		};
		DefaultUnaryBinaryReactionNetwork net = new DefaultUnaryBinaryReactionNetwork(productionStochiometries[0].length, productionStochiometries.length);
		net.setStoichiometries(productionStochiometries, consumptionStochiometries);
		net.setRateParameters(rateParameters);
		String[] speciesNames = { "S" };

		this.net = net;
		this.deterministicReactions = deterministicReactions;
		this.N = N;
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
