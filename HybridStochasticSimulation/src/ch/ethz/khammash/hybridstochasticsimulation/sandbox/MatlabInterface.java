package ch.ethz.khammash.hybridstochasticsimulation.sandbox;

import java.lang.reflect.Array;

import ch.ethz.khammash.hybridstochasticsimulation.models.MSHybridReactionNetworkModel;
import ch.ethz.khammash.hybridstochasticsimulation.networks.DefaultUnaryBinaryReactionNetwork;
import ch.ethz.khammash.hybridstochasticsimulation.networks.MSHybridReactionNetwork;


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

//	public PDMPModelAdapter createPDMPModel(FirstOrderDifferentialEquations deterministicModel,
//			StochasticReactionNetworkModel stochasticModel) {
//		return new PDMPModelAdapter(new HybridModelAdapter(deterministicModel, stochasticModel));
//	}
//
//	public PDMPSimulator createPDMPModelSimulator() {
//		return new PDMPSimulator();
//	}
//
//	public ContinuousPDMPTrajectory createPDMPContinuousOutputModel() {
//		return new ContinuousPDMPTrajectory();
//	}

	public Object createDoubleArray(int length) {
		return Array.newInstance(double[].class, length);
	}

	public Object createDoubleArray(int... dimensions) {
		return Array.newInstance(double[].class, dimensions);
	}

}
