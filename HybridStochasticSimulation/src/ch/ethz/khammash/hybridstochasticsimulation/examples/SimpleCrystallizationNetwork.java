package ch.ethz.khammash.hybridstochasticsimulation.examples;

import org.apache.commons.math3.random.MersenneTwister;
import org.apache.commons.math3.random.RandomDataGenerator;

import ch.ethz.khammash.hybridstochasticsimulation.networks.ReactionNetwork;

public class SimpleCrystallizationNetwork extends ExampleNetwork {

	public SimpleCrystallizationNetwork() {
		int[] continuousSpecies = { };
		double N = 1e6;
		double deltaR = 0.5;
		double deltaS = 0.5;
		double epsilon = 0.5;
		double gamma = 0;
		double[] alpha = { 1, 1, 0, 0 };
		double[] beta = { -1, -1 };
		double t0 = 0.0;
		double t1 = 100.0;
		double[] x0 = { 1e6, 0, 10, 0 };
		double[] plotScales = { 1e-5, 1e-5, 1, 1 };
		ReactionNetwork net = new ReactionNetwork(4, 2);
		net.setStochiometry(0, 0, 0, 2);
		net.setStochiometry(1, 0, 1, 0);
		net.setStochiometry(2, 0, 0, 0);
		net.setStochiometry(3, 0, 0, 0);
		net.setStochiometry(0, 1, 0, 1);
		net.setStochiometry(1, 1, 0, 0);
		net.setStochiometry(2, 1, 0, 1);
		net.setStochiometry(3, 1, 1, 0);
		net.setRateParameter(0, 1e-7);
		net.setRateParameter(1, 1e-7);
		String[] speciesNames = { "A", "B", "C", "D" };

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
