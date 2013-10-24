package ch.ethz.khammash.hybridstochasticsimulation.simulators;

public class PDMPSimulationInformation {

	private long integrationSteps = 0;
	private boolean integrating = false;
	private long totalReactionCount = 0;
	private long[] reactionCounts;

	public PDMPSimulationInformation(int numOfReactions) {
		reactionCounts = new long[numOfReactions];
	}

	public long getIntegrationSteps() {
		return integrationSteps;
	}

	public long getTotalReactionCount() {
		return totalReactionCount;
	}

	public double getRelativeReactionCount(int reaction) {
		return getReactionCount(reaction) / (double)getTotalReactionCount();
	}

	public double[] getRelativeReactionCounts() {
		double[] arr = new double[reactionCounts.length];
		for (int i=0; i < reactionCounts.length; i++)
			arr[i] = getRelativeReactionCount(i);
		return arr;
	}

	public long getReactionCount(int reaction) {
		return reactionCounts[reaction];
	}

	public long[] getReactionCounts() {
		return reactionCounts.clone();
	}

	public boolean isIntegrating() {
		return integrating;
	}

	public double[] computeInformationState() {
		double[] state = new double[3 + 2 * reactionCounts.length];
		state[0] = isIntegrating() ? 1.0 : -1.0;
		state[1] = getIntegrationSteps();
		state[2] = getTotalReactionCount();
		for (int i=0; i < reactionCounts.length; i++)
			state[3 + i] = getReactionCount(i);
		for (int i=0; i < reactionCounts.length; i++)
			state[3 + reactionCounts.length + i] = getRelativeReactionCount(i);
		return state;
	}

	public void setIntegrationOn() {
		integrating = true;
	}

	public void setIntegrationOff() {
		integrating = false;
	}

	public void increaseIntegrationSteps() {
		integrationSteps++;
	}

	public void increaseReactionCount(int reaction) {
		totalReactionCount++;
		reactionCounts[reaction]++;
	}

}
