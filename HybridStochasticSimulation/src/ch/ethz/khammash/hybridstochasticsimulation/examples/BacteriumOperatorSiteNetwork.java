package ch.ethz.khammash.hybridstochasticsimulation.examples;

import org.apache.commons.math3.random.MersenneTwister;
import org.apache.commons.math3.random.RandomDataGenerator;

import ch.ethz.khammash.hybridstochasticsimulation.models.ReactionNetwork;



public class BacteriumOperatorSiteNetwork extends ExampleNetwork {

	public BacteriumOperatorSiteNetwork() {
		int[] continuousSpecies = { 2 };
		double N = 100;
		double deltaR = 0.5;
		double deltaS = 0.5;
		double epsilon = 0.1;
		double gamma = 0;
		double[] alpha = { 0, 0, 1 };
		double[] beta = { 0, 0, 1, 0 };
		double t0 = 0.0;
		double t1 = 40 * 60 * 60;
		double[] x0 = { 1, 0, 0 };
		double[] plotScales = { 100, 100, 1 };
		int[][] productionStochiometries = {
				{ 0, 1, 0 },
				{ 1, 0, 0 },
				{ 0, 1, 1 },
				{ 0, 0, 0 },
		};
		int[][] consumptionStochiometries = {
				{ 1, 0, 0 },
				{ 0, 1, 0 },
				{ 0, 1, 0 },
				{ 0, 0, 1 },
		};
		double[] rateParameters = {
				2.0e-3,
				2.0e-3,
				1,
				4.8e-4,
		};
		ReactionNetwork net = new ReactionNetwork(productionStochiometries[0].length, productionStochiometries.length);
		net.setStochiometries(productionStochiometries, consumptionStochiometries);
		net.setRateParameters(rateParameters);
		String[] speciesNames = { "G", "G*", "P" };

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
