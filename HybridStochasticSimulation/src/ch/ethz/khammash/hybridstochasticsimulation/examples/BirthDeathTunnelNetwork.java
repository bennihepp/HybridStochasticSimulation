package ch.ethz.khammash.hybridstochasticsimulation.examples;

import org.apache.commons.math3.random.MersenneTwister;
import org.apache.commons.math3.random.RandomDataGenerator;

import ch.ethz.khammash.hybridstochasticsimulation.networks.DefaultUnaryBinaryReactionNetwork;



public class BirthDeathTunnelNetwork extends ExampleNetwork {

	public BirthDeathTunnelNetwork() {
		// See Kang & Kurtz 2011 (3.2 Collective species balance)
		int[] continuousSpecies = { };
		double N = 1e6;
		double deltaR = 0.5;
		double deltaS = 0.5;
		double epsilon = 0.5;
		double gamma = 0;
		double[] alpha = { 0, 0 };
		double[] beta = { 0, 0, 0, 0 };
		double t0 = 0.0;
		double t1 = 100.0;
		double[] x0 = { 0, 0 };
		double[] plotScales = { 1, 1 };
		int[][] productionStochiometries = {
				{ 1, 0 },
				{ 0, 1 },
				{ 1, 0 },
				{ 0, 0 },
		};
		int[][] consumptionStochiometries = {
				{ 0, 0 },
				{ 1, 0 },
				{ 0, 1 },
				{ 0, 1 },
		};
//		double[] rateParameters = {
//				1e-1,
//				1.0,
//				1.0,
//				1e-2,
//		};
		double[] rateParameters = {
				1,
				10.0,
				10.0,
				1e-1,
		};
		DefaultUnaryBinaryReactionNetwork net = new DefaultUnaryBinaryReactionNetwork(productionStochiometries[0].length, productionStochiometries.length);
		net.setStochiometries(productionStochiometries, consumptionStochiometries);
		net.setRateParameters(rateParameters);
		String[] speciesNames = { "S1", "S2" };

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