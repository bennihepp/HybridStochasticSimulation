package ch.ethz.bhepp.hybridstochasticsimulation.examples;

import org.apache.commons.math3.random.MersenneTwister;
import org.apache.commons.math3.random.RandomDataGenerator;

import ch.ethz.bhepp.hybridstochasticsimulation.networks.DefaultUnaryBinaryReactionNetwork;



// Haploinsufficiency Example from:
// Crudu, Alina, Arnaud Debussche, and Ovidiu Radulescu. "Hybrid stochastic simplifications for multiscale gene networks." BMC systems biology 3.1 (2009): 89.
//
// species:
//  S0: G
//  S1: G*
//  S2: P
// reactions:
//  R0: G -> G*       [2e-4]
//  R1: G* -> G       [2e-4]
//  R2: G* -> G* + P  [10]
//  R3: P -> ~        [4.8e-4]

public class HaploinsufficiencyNetwork extends SimulationConfiguration {

	public HaploinsufficiencyNetwork() {
		int[] deterministicReactions = { 2 };
		double N = 100;
		double gamma = 0;
		double[] alpha = { 0, 0, 1 };
		double[] beta = { 0, 0, 1, 0 };
		double t0 = 0.0;
		double t1 = 40 * 60 * 60;
		double[] x0 = { 1, 0, 0 };
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
//		double[] rateParameters = {
//				2.0e-3,
//				2.0e-3,
//				1,
//				4.8e-4,
//		};
		double[] rateParameters = {
				2.0e-4,
				2.0e-4,
				10,
				4.8e-4,
		};
		double[] plotScales = { (rateParameters[2] / rateParameters[3]) / 5, (rateParameters[2] / rateParameters[3]) / 5, 1  };
		DefaultUnaryBinaryReactionNetwork net = new DefaultUnaryBinaryReactionNetwork(productionStochiometries[0].length, productionStochiometries.length);
		net.setStoichiometries(productionStochiometries, consumptionStochiometries);
		net.setRateParameters(rateParameters);
		String[] speciesNames = { "G", "G*", "P" };

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
