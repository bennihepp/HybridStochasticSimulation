package ch.ethz.khammash.hybridstochasticsimulation.batch.guiceproviders;

import org.apache.commons.configuration.HierarchicalConfiguration;

import ch.ethz.khammash.hybridstochasticsimulation.trajectories.ArrayFiniteTrajectoryRecorder;
import ch.ethz.khammash.hybridstochasticsimulation.trajectories.FiniteTrajectoryRecorder;

import com.google.inject.Inject;

public class FiniteTrajectoryRecorderProvider extends AbstractObjProvider<FiniteTrajectoryRecorder> {

	@Inject
	public FiniteTrajectoryRecorderProvider(HierarchicalConfiguration config) {
		super(config, "SimulationParameters");
	}

	@Override
	public FiniteTrajectoryRecorder get() {
		int numOfTimePoints = config().getInt("numOfTimePoints");
		return new ArrayFiniteTrajectoryRecorder(numOfTimePoints);
	}

}
