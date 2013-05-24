package ch.ethz.khammash.hybridstochasticsimulation;

import java.lang.reflect.Array;

import org.apache.commons.math3.ode.FirstOrderDifferentialEquations;

import ch.ethz.khammash.hybridstochasticsimulation.models.MSHybridReactionNetworkModel;
import ch.ethz.khammash.hybridstochasticsimulation.models.PDMPModelAdapter;
import ch.ethz.khammash.hybridstochasticsimulation.models.PDMPModelTrajectory;
import ch.ethz.khammash.hybridstochasticsimulation.models.ReactionNetworkModel;
import ch.ethz.khammash.hybridstochasticsimulation.networks.MSHybridReactionNetwork;
import ch.ethz.khammash.hybridstochasticsimulation.networks.ReactionNetwork;
import ch.ethz.khammash.hybridstochasticsimulation.simulators.PDMPModelSimulator;


public class MatlabInterface {

	public ReactionNetwork createReactionNetwork(int numOfSpecies,
			int numOfReactions) {
		return new ReactionNetwork(numOfSpecies, numOfReactions);
	}

	public MSHybridReactionNetwork createHybridReactionNetwork(
			ReactionNetwork net, double N, double deltaR, double deltaS, double gamma, double[] alpha, double[] beta) {
		return new MSHybridReactionNetwork(net, N, deltaR, deltaS, gamma, alpha, beta);
	}

	public MSHybridReactionNetworkModel createHybridReactionNetworkModel(
			MSHybridReactionNetwork net) {
		return new MSHybridReactionNetworkModel(net);
	}

	public PDMPModelAdapter createPDMPModel(FirstOrderDifferentialEquations baseODE,
			ReactionNetworkModel reactionModel) {
		return new PDMPModelAdapter(baseODE, reactionModel);
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
