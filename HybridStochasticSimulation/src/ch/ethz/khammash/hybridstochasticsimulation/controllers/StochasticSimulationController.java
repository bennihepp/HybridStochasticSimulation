package ch.ethz.khammash.hybridstochasticsimulation.controllers;

import java.util.concurrent.ExecutorService;

import org.apache.commons.math3.random.RandomDataGenerator;

import ch.ethz.khammash.hybridstochasticsimulation.models.StochasticReactionNetworkModel;
import ch.ethz.khammash.hybridstochasticsimulation.simulators.Simulator;
import ch.ethz.khammash.hybridstochasticsimulation.simulators.StochasticSimulator;
import ch.ethz.khammash.hybridstochasticsimulation.trajectories.TrajectoryRecorder;


public class StochasticSimulationController<T extends StochasticReactionNetworkModel>
		extends AbstractSimulationController<T, TrajectoryRecorder<T>> {

    public StochasticSimulationController() {
    	super();
		construct();
    }

	public StochasticSimulationController(int numOfThreads) {
		super(numOfThreads);
		construct();
    }

	public StochasticSimulationController(ExecutorService executor) {
		super(executor);
		construct();
    }

	final private void construct() {
		setSimulatorFactory(
				new SimulatorFactory<Simulator<T, TrajectoryRecorder<T>>>() {

			@Override
			public Simulator<T, TrajectoryRecorder<T>> createSimulator(
					RandomDataGenerator rdg) {
				return new StochasticSimulator<T>(rdg);
			}

		});
//		setTrajectoryRecorderFactory(new TrajectoryRecorderFactory<TrajectoryRecorder<StochasticReactionNetworkModel>>() {
//
//			@Override
//			public TrajectoryRecorder<StochasticReactionNetworkModel> createTrajectoryRecorder() {
//				return new 
//			}
//
//		})
	}

}
