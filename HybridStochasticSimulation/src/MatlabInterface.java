import java.lang.reflect.Array;

import org.apache.commons.math3.ode.FirstOrderDifferentialEquations;


public class MatlabInterface {

	public ReactionNetwork createReactionNetwork(int numOfSpecies, int numOfReactions) {
		return new ReactionNetwork(numOfSpecies, numOfReactions);
	}

	public HybridReactionNetwork createHybridReactionNetwork(
			ReactionNetwork net, double N, int gamma, int[] alpha, int[] beta) {
		return new HybridReactionNetwork(net, N, gamma, alpha, beta);
	}

	public HybridReactionNetworkModel createHybridReactionNetworkModel(HybridReactionNetwork net) {
		return new HybridReactionNetworkModel(net);
	}

	public PDMPModel createPDMPModel(FirstOrderDifferentialEquations baseODE, ReactionNetworkModel reactionModel) {
		return new PDMPModel(baseODE, reactionModel);
	}

	public PDMPModelSimulator createPDMPModelSimulator() {
		return new PDMPModelSimulator();
	}

	public PDMPModelTrajectory createPDMPContinuousOutputModel() {
		return new PDMPModelTrajectory();
	}

	public Object createDoubleArray(int length) {
		return Array.newInstance(double[].class, length);
	}

	public Object createDoubleArray(int... dimensions) {
		return Array.newInstance(double[].class, dimensions);
	}

}
