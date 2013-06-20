package ch.ethz.khammash.hybridstochasticsimulation.examples;

import org.apache.commons.math3.random.MersenneTwister;
import org.apache.commons.math3.random.RandomDataGenerator;

import ch.ethz.khammash.hybridstochasticsimulation.networks.DefaultUnaryBinaryReactionNetwork;



public class LambdaPhageToggleSwitchNetwork extends ExampleConfiguration {

	public LambdaPhageToggleSwitchNetwork() {
		// See Crudu et. al 2009
		// species:
		//  D
		//  DcI2
		//  DcI2*
		//  DcI2cI2
		//  cI
		//  cI2
		// reactions:
		//  2cI -> cI2 (k1)
		//  cI2 -> 2cI (km1)
		//  DcI2 -> DcI2 + ncI (k5)
		//  cI -> - (k6)
		//  D + cI2 -> DcI2 (k2)
		//  DcI2 -> D + cI2 (km2)
		//  D + cI2 -> DcI2* (k3)
		//  DcI2* -> D + cI2 (km3)
		//  DcI2 + cI2 -> DcI2cI2 (k4)
		//  DcI2cI2 -> DcI2 + cI2 (km4)
		int n = 10;
		int[] deterministicReactions = {};
		double N = 100;
		double gamma = 0;
		double[] alpha = { 0, 0, 0 };
		double[] beta = { 0, 0, 0, 0 };
		double t0 = 0.0;
		double t1 = 2e4;
		double[] x0 = { 1, 0, 0, 0, 50, 50 };
		double[] plotScales = { 1, 1, 1, 1, 1, 1 };
		int[][] productionStochiometries = {
				{ 0, 0, 0, 0, 0, 1 }, // 2cI -> cI2 (k1)
				{ 0, 0, 0, 0, 2, 0 }, // cI2 -> 2cI (km1)
				{ 0, 1, 0, 0, n, 0 }, // DcI2 -> DcI2 + ncI (k5)
				{ 0, 0, 0, 0, 0, 0 }, // cI -> - (k6)
				{ 0, 1, 0, 0, 0, 0 }, // D + cI2 -> DcI2 (k2)
				{ 1, 0, 0, 0, 0, 1 }, // DcI2 -> D + cI2 (km2)
				{ 0, 0, 1, 0, 0, 0 }, // D + cI2 -> DcI2* (k3)
				{ 1, 0, 0, 0, 0, 1 }, // DcI2* -> D + cI2 (km3)
				{ 0, 0, 0, 1, 0, 0 }, // DcI2 + cI2 -> DcI2cI2 (k4)
				{ 0, 1, 0, 0, 0, 1 }, // DcI2cI2 -> DcI2 + cI2 (km4)
		};
		int[][] consumptionStochiometries = {
				{ 0, 0, 0, 0, 2, 0 }, // 2cI -> cI2 (k1)
				{ 0, 0, 0, 0, 0, 1 }, // cI2 -> 2cI (km1)
				{ 0, 1, 0, 0, 0, 0 }, // DcI2 -> DcI2 + ncI (k5)
				{ 0, 0, 0, 0, 1, 0 }, // cI -> - (k6)
				{ 1, 0, 0, 0, 0, 1 }, // D + cI2 -> DcI2 (k2)
				{ 0, 1, 0, 0, 0, 0 }, // DcI2 -> D + cI2 (km2)
				{ 1, 0, 0, 0, 0, 1 }, // D + cI2 -> DcI2* (k3)
				{ 0, 0, 1, 0, 0, 0 }, // DcI2* -> D + cI2 (km3)
				{ 0, 1, 0, 0, 0, 1 }, // DcI2 + cI2 -> DcI2cI2 (k4)
				{ 0, 0, 0, 1, 0, 0 }, // DcI2cI2 -> DcI2 + cI2 (km4)
		};
		double k1 =  1.0;
		double km1 = 1.0;
		double k5 =  1.0;
		double k6 =  1.0;
		double k2 =  1.0;
		double km2 = 1.0;
		double k3 =  1.0;
		double km3 = 1.0;
		double k4 =  1.0;
		double km4 = 1.0;
		double[] rateParameters = {
				k1,
				km1,
				k5,
				k6,
				k2,
				km2,
				k3,
				km3,
				k4,
				km4,
		};
		DefaultUnaryBinaryReactionNetwork net = new DefaultUnaryBinaryReactionNetwork(productionStochiometries[0].length, productionStochiometries.length);
		net.setStochiometries(productionStochiometries, consumptionStochiometries);
		net.setRateParameters(rateParameters);
		String[] speciesNames = { "D", "DcI2", "DcI2*", "DcI2cI2", "cI", "cI2" };

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
