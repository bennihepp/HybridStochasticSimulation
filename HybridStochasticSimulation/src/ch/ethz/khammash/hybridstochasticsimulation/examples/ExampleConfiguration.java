package ch.ethz.khammash.hybridstochasticsimulation.examples;


import org.apache.commons.math3.random.RandomDataGenerator;
import org.apache.commons.math3.random.RandomGenerator;

import ch.ethz.khammash.hybridstochasticsimulation.networks.DefaultUnaryBinaryReactionNetwork;


public class ExampleConfiguration {
	public RandomGenerator rng;
	public RandomDataGenerator rdg;
	public DefaultUnaryBinaryReactionNetwork net;
	public double N;
	public double gamma;
	public double[] alpha;
	public double[] beta;
	public double t0;
	public double t1;
	public double[] x0;
	public double[] plotScales;
	public String[] speciesNames;
	public int[] deterministicReactions = { };
	public double epsilon = 0.1;
	public double delta = 0.5;
	public double xi = 0.5;
	public double theta = 1000;
	// How many reactions should occur before the copy number bounds are checked manually
	public int zeta = 100;
	public double tolerance = 1e-6;
}
