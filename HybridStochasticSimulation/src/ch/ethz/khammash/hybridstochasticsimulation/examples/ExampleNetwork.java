package ch.ethz.khammash.hybridstochasticsimulation.examples;


import org.apache.commons.math3.random.RandomDataGenerator;
import org.apache.commons.math3.random.RandomGenerator;

import ch.ethz.khammash.hybridstochasticsimulation.networks.DefaultUnaryBinaryReactionNetwork;


public class ExampleNetwork {
	public RandomGenerator rng;
	public RandomDataGenerator rdg;
	public DefaultUnaryBinaryReactionNetwork net;
	public int[] continuousSpecies;
	public double N;
	public double deltaR;
	public double deltaS;
	public double epsilon;
	public double gamma;
	public double[] alpha;
	public double[] beta;
	public double t0;
	public double t1;
	public double[] x0;
	public double[] plotScales;
	public String[] speciesNames;
}