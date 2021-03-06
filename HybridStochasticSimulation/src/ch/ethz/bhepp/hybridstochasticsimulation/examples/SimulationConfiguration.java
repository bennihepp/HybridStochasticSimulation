package ch.ethz.bhepp.hybridstochasticsimulation.examples;import java.util.Arrays;

import org.apache.commons.math3.random.RandomDataGenerator;
import org.apache.commons.math3.random.RandomGenerator;

import ch.ethz.bhepp.hybridstochasticsimulation.networks.MassActionReactionNetwork;


public class SimulationConfiguration {

	public RandomGenerator rng;
	public RandomDataGenerator rdg;
	public MassActionReactionNetwork net;
	public double N;
	public double gamma;
	public double[] alpha;
	public double[] beta;
	public double t0;
	public double tf;
	public double[] x0;
	public double[] plotScales;
	public String[] speciesNames;
	public int[] importantSpecies = {};
	public int[] deterministicReactions = {};
	public double eta = 0.1;
	public double delta = 0.5;
	public double mu = 0.5;
	public double theta = 1;
	// How many reactions should occur before the copy number bounds are checked manually
	public int zeta = 100;
	public double tolerance = 1e-6;

	public void reset() {
		if (net != null) {
			x0 = new double[net.getNumberOfSpecies()];
			plotScales = new double[net.getNumberOfSpecies()];
			Arrays.fill(plotScales, 1.0);
			speciesNames = new String[net.getNumberOfSpecies()];
			for (int s=0; s < net.getNumberOfSpecies(); s++)
				speciesNames[s] = "S" + s;
			importantSpecies = new int[0];
			deterministicReactions = new int[0];
			alpha = new double[net.getNumberOfSpecies()];
			beta = new double[net.getNumberOfReactions()];
		}
	}
}
