package ch.ethz.khammash.hybridstochasticsimulation.simulators;

import static com.google.common.base.Preconditions.checkArgument;

import java.util.LinkedList;
import java.util.List;

import org.apache.commons.math3.random.RandomDataGenerator;

import ch.ethz.khammash.hybridstochasticsimulation.Utilities;
import ch.ethz.khammash.hybridstochasticsimulation.math.RandomDataUtilities;
import ch.ethz.khammash.hybridstochasticsimulation.models.StochasticReactionNetworkModel;
import ch.ethz.khammash.hybridstochasticsimulation.trajectories.TrajectoryRecorder;


public class StochasticSimulator
		extends AbstractSimulator<StochasticReactionNetworkModel> {

	private RandomDataGenerator rdg;
	private List<TrajectoryRecorder> trajectoryRecorders;

	public StochasticSimulator() {
		this(null);
	}

	public StochasticSimulator(RandomDataGenerator rdg) {
		if (rdg == null)
			rdg = new RandomDataGenerator();
		this.rdg = rdg;
		trajectoryRecorders = new LinkedList<TrajectoryRecorder>();
	}

	public double simulate(StochasticReactionNetworkModel model, double t0, double[] x0, double t1, double[] x1) {
		checkArgument(x0.length == x1.length, "Expected x0.length == x1.length");
		checkArgument(x0.length == model.getNumberOfSpecies(), "Expected x0.length == model.getNumberOfSpecies()");

		final int STOCHASTIC_RECORD_INTERVAL = 1;

		double[] x = new double[x0.length];
		for (int i=0; i < x0.length; i++)
			x[i] = x0[i];
		double t = t0;
		double[] propVec = new double[model.getNumberOfReactions()];
    	for (TrajectoryRecorder handler : trajectoryRecorders) {
    		handler.beginRecording(t0, x0, t1);
    	}
    	long reactionCounter = 0;
    	double[] reactionCounterArray = new double[model.getNumberOfReactions()];
    	double msgDt = (t1 - t0) / 20.0;
    	double nextMsgT = t0 + msgDt;
		final long startTime = System.currentTimeMillis();

		int j = 0;

		while (true) {

			if (showProgress)
				while (t > nextMsgT) {
					System.out.println("Progress: " + (100 * (nextMsgT - t0)/(t1 - t0)) + "%");
					nextMsgT += msgDt;
				}

	        double propSum = model.computePropensitiesAndSum(t, x, propVec);

	        // Sanity check for negative propensities or no more reactions occuring (zero propensities)
	        if (propSum < 0)
	        	throw new RuntimeException("Negative propensities are not allowed to occur!");
	        else if (propSum == 0)
	        	break;

	        double tau = rdg.nextExponential(1 / propSum);
	        t = t + tau;

	        // Stop if we reached the end-timepoint
	        if (t >= t1)
	        	break;

	        // Determine which reaction fired and update state
	        int reaction = RandomDataUtilities.sampleFromProbabilityMassFunction(rdg, propVec, propSum);
    		model.changeState(reaction, t, x);

        	reactionCounter++;
        	reactionCounterArray[reaction]++;

        	if (j > STOCHASTIC_RECORD_INTERVAL) {
	        	for (TrajectoryRecorder handler : trajectoryRecorders)
	        		handler.record(t, x);
        		j = 0;
        	}
        	j++;

		}

		final long endTime = System.currentTimeMillis();
		if (printMessages) {
			System.out.println("Execution time: " + (endTime - startTime));
			System.out.println("Total of " + reactionCounter + " reactions performed");
			Utilities.printArray("Total reaction counts", reactionCounterArray);
			for (int r=0; r < reactionCounterArray.length; r++)
				reactionCounterArray[r] = 100.0 * reactionCounterArray[r] / reactionCounter;
			Utilities.printArray("Relative reaction counts", reactionCounterArray, "%.2f");
		}
		for (int i=0; i < x1.length; i++)
			x1[i] = x[i];
    	for (TrajectoryRecorder handler : trajectoryRecorders)
    		handler.endRecording(x1);
		return t;
	}

	@Override
	public void addTrajectoryRecorder(TrajectoryRecorder tr) {
		trajectoryRecorders.add(tr);
	}

	@Override
	public void removeTrajectoryRecorder(TrajectoryRecorder tr) {
		trajectoryRecorders.remove(tr);
	}

	@Override
	public void clearTrajectoryRecorders() {
		trajectoryRecorders.clear();
	}

}
