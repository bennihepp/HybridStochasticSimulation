package ch.ethz.khammash.hybridstochasticsimulation.sandbox;

import java.lang.reflect.Array;


import org.apache.commons.math3.ode.FirstOrderDifferentialEquations;

import ch.ethz.khammash.hybridstochasticsimulation.models.HybridModelAdapter;
import ch.ethz.khammash.hybridstochasticsimulation.models.MSHybridReactionNetworkModel;
import ch.ethz.khammash.hybridstochasticsimulation.models.PDMPModelAdapter;
import ch.ethz.khammash.hybridstochasticsimulation.models.StochasticReactionNetworkModel;
import ch.ethz.khammash.hybridstochasticsimulation.networks.MSHybridReactionNetwork;
import ch.ethz.khammash.hybridstochasticsimulation.networks.DefaultUnaryBinaryReactionNetwork;
import ch.ethz.khammash.hybridstochasticsimulation.simulators.PDMPSimulator;
import ch.ethz.khammash.hybridstochasticsimulation.trajectories.PDMPTrajectory;


public class MatlabInterface {

	public DefaultUnaryBinaryReactionNetwork createReactionNetwork(int numOfSpecies,
			int numOfReactions) {
		return new DefaultUnaryBinaryReactionNetwork(numOfSpecies, numOfReactions);
	}

	public MSHybridReactionNetwork createHybridReactionNetwork(
			DefaultUnaryBinaryReactionNetwork net, double N, double deltaR, double deltaS, double epsilon, double gamma, double[] alpha, double[] beta) {
		return new MSHybridReactionNetwork(net, N, deltaR, deltaS, epsilon, gamma, alpha, beta);
	}

	public MSHybridReactionNetworkModel createHybridReactionNetworkModel(
			MSHybridReactionNetwork net) {
		return new MSHybridReactionNetworkModel(net);
	}

	public PDMPModelAdapter createPDMPModel(FirstOrderDifferentialEquations deterministicModel,
			StochasticReactionNetworkModel stochasticModel) {
		return new PDMPModelAdapter(new HybridModelAdapter(deterministicModel, stochasticModel));
	}

	public PDMPSimulator createPDMPModelSimulator() {
		return new PDMPSimulator();
	}

	public PDMPTrajectory createPDMPContinuousOutputModel() {
		return new PDMPTrajectory();
	}

	public Object createDoubleArray(int length) {
		return Array.newInstance(double[].class, length);
	}

	public Object createDoubleArray(int... dimensions) {
		return Array.newInstance(double[].class, dimensions);
	}

}
