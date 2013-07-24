package ch.ethz.khammash.hybridstochasticsimulation.batch;

import ch.ethz.khammash.hybridstochasticsimulation.trajectories.FiniteTrajectory;
import de.root1.simon.annotation.SimonRemote;

@SimonRemote(value = {ClientCallbackInterface.class}) 
public class ClientCallbackImpl implements ClientCallbackInterface {

	@Override
	public void callback(FiniteTrajectory tr) {
		System.out.println("Trajectory: " + tr);
	}

}
