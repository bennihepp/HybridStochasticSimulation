package ch.ethz.khammash.hybridstochasticsimulation.batch.providers;

import org.apache.commons.configuration.HierarchicalConfiguration;

import ch.ethz.khammash.hybridstochasticsimulation.factories.FiniteTrajectoryRecorderFactory;
import ch.ethz.khammash.hybridstochasticsimulation.trajectories.ArrayFiniteTrajectoryRecorder;
import ch.ethz.khammash.hybridstochasticsimulation.trajectories.FiniteTrajectoryRecorder;

import com.google.inject.Inject;

public class FiniteTrajectoryRecorderFactoryProvider extends AbstractProvider<FiniteTrajectoryRecorderFactory> {

	@Inject
	public FiniteTrajectoryRecorderFactoryProvider(HierarchicalConfiguration config) {
		super(config, "SimulationParameters");
	}

	@Override
	public FiniteTrajectoryRecorderFactory get() {
		final int numOfTimePoints = config().getInt("numOfTimePoints");
		return new FiniteTrajectoryRecorderFactory() {
			
			@Override
			public FiniteTrajectoryRecorder createTrajectoryRecorder() {
				return new ArrayFiniteTrajectoryRecorder(numOfTimePoints);
			}

		};
	}

}
