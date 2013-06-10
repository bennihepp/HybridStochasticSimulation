package ch.ethz.khammash.hybridstochasticsimulation.examples;


import org.apache.commons.math3.random.MersenneTwister;
import org.apache.commons.math3.random.RandomDataGenerator;

import ch.ethz.khammash.hybridstochasticsimulation.models.ReactionNetwork;


public class StochasticFocusingNetwork extends ExampleNetwork {

	public StochasticFocusingNetwork() {
		// See Paulsson et al. PNAS 2000
		// species:
		//  I
		//  P
		//  S
		// reactions:
		//  - -> I (k)
		//  I -> - (k*S)
		//  I -> P (kp)
		//  P -> - (1)
		//  - -> S (ks)
		//  S -> - (kd)
		int[] continuousSpecies = { };
		double N = 1e4;
		double deltaR = 0.5;
		double deltaS = 0.5;
		double epsilon = 0.5;
		double gamma = 0;
		double[] alpha = { 0, 0, 0 };
		double[] beta = { 0, 0, 0, 0, 0, 0 };
		double t0 = 0.0;
		double t1 = 100.0;
		double[] x0 = { 0, 0, 0 };
		double[] plotScales = { 1, 1, 1 };
		int[][] productionStochiometries = {
				{ 1, 0, 0 }, //  - -> I (k)
				{ 0, 0, 1 }, //  I -> - (ka*S)
				{ 0, 1, 0 }, //  I -> P (kp)
				{ 0, 0, 0 }, //  P -> - (1)
				{ 0, 0, 1 }, //  - -> S (ks)
				{ 0, 0, 0 }, //  S -> - (kd)
		};
		int[][] consumptionStochiometries = {
				{ 0, 0, 0 }, //  - -> I (k)
				{ 1, 0, 1 }, //  I -> - (ka*S)
				{ 1, 0, 0 }, //  I -> P (kp)
				{ 0, 1, 0 }, //  P -> - (1)
				{ 0, 0, 0 }, //  - -> S (ks)
				{ 0, 0, 1 }, //  S -> - (kd)
		};
//		double k = 1e5;
//		double ka = 100.0;
//		double kp = 100;
//		double kd = 1.0;
//		double ks;
////		ks = 10.0 * kd;
//		ks = 5.0 * kd;
		double k = 1e6;
		double ka = 100.0;
		double kp = 100;
		double kd = 10.0;
		double ks;
//		ks = 10.0 * kd;
		ks = 5.0 * kd;
		double[] rateParameters = {
				k,
				ka,
				kp,
				1.0,
				ks,
				kd,
		};
		ReactionNetwork net = new ReactionNetwork(productionStochiometries[0].length, productionStochiometries.length);
		net.setStochiometries(productionStochiometries, consumptionStochiometries);
		net.setRateParameters(rateParameters);
		String[] speciesNames = { "I", "P", "S" };

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
