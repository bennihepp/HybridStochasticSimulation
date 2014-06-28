package ch.ethz.bhepp.hybridstochasticsimulation.grid.akka;

import java.util.HashMap;
import java.util.Map;

import akka.actor.ActorRef;
import akka.actor.UntypedActor;
import ch.ethz.bhepp.hybridstochasticsimulation.grid.akka.AkkaWorker.WorkerMessage;

public class AkkaController extends UntypedActor {

	public static enum ControllerMessage {
		RUN_SIMULATION, SHUTDOWN,
	}

	private Map<ActorRef, Integer> workers;

	public AkkaController() {
		this.workers = new HashMap<>();
	}

	@Override
	public void onReceive(Object msg) throws Exception {
		if (msg instanceof WorkerMessage) {
			if (msg == WorkerMessage.IDLE) {
				workers.put(sender(), workers.size());
				System.out.println("Worker is waiting");
				sender().tell(ControllerMessage.RUN_SIMULATION, self());
			}
		}
	}

}
