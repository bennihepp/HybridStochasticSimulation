package ch.ethz.khammash.hybridstochasticsimulation.examples;

import org.apache.commons.math3.random.MersenneTwister;
import org.apache.commons.math3.random.RandomDataGenerator;

import ch.ethz.khammash.hybridstochasticsimulation.networks.DefaultUnaryBinaryReactionNetwork;



public class RegulatedTranscriptionNetwork extends ExampleConfiguration {

	public RegulatedTranscriptionNetwork() {
		int[] deterministicReactions = { };
		double N = 100;
		double gamma = 0;
		double[] alpha = { 1, 1, 0, 0, 0, 0 };
		double[] beta = { -1, -2, -1, -1, -1, 0, -3, -2, -1, 0 };
		double t0 = 0.0;
		double t1 = 100.0;
		double[] x0 = { 2, 6, 0, 0, 2, 0 };
		double[] plotScales = { 1, 1, 1, 1, 1, 1 };
		int[][] productionStochiometries = {
				{ 1, 0, 1, 0, 0, 0 },
				{ 0, 0, 0, 0, 0, 0 },
				{ 0, 0, 1, 0, 1, 0 },
				{ 0, 0, 0, 0, 0, 0 },
				{ 0, 0, 0, 0, 1, 0 },
				{ 0, 1, 0, 1, 0, 0 },
				{ 0, 0, 0, 0, 0, 1 },
				{ 0, 1, 0, 0, 1, 0 },
				{ 0, 1, 0, 0, 0, 0 },
				{ 2, 0, 0, 0, 0, 0 }
		};
		int[][] consumptionStochiometries = {
				{ 0, 0, 1, 0, 0, 0 },
				{ 1, 0, 0, 0, 0, 0 },
				{ 0, 0, 0, 0, 1, 0 },
				{ 0, 0, 1, 0, 0, 0 },
				{ 0, 1, 0, 1, 0, 0 },
				{ 0, 0, 0, 0, 1, 0 },
				{ 0, 1, 0, 0, 1, 0 },
				{ 0, 0, 0, 0, 0, 1 },
				{ 2, 0, 0, 0, 0, 0 },
				{ 0, 1, 0, 0, 0, 0 }
		};
		double[] rateParameters = {
				4.30e-2,
				7.00e-4,
				7.15e-2,
				3.90e-3,
				1.99e-2,
				4.79e-1,
				1.99e-4,
				8.77e-12,
				8.30e-2,
				5.00e-1
		};
		DefaultUnaryBinaryReactionNetwork net = new DefaultUnaryBinaryReactionNetwork(productionStochiometries[0].length, productionStochiometries.length);
		net.setStochiometries(productionStochiometries, consumptionStochiometries);
		net.setRateParameters(rateParameters);
		String[] speciesNames = { "M", "D", "RNA", "DNA", "DNAD", "DNA2D" };

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
