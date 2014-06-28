package ch.ethz.bhepp.hybridstochasticsimulation.grid.akka;

import akka.actor.ActorIdentity;
import akka.actor.ActorRef;
import akka.actor.Identify;
import akka.actor.ReceiveTimeout;
import akka.actor.UntypedActor;
import ch.ethz.bhepp.hybridstochasticsimulation.grid.akka.AkkaController.ControllerMessage;

public class AkkaWorker extends UntypedActor {

	public static enum WorkerMessage {
		IDLE,
	}

	private ActorRef master;
	private String path;

	public AkkaWorker(String path) {
		this.path = path;
		sendIdentifyRequest();
	}

	private void sendIdentifyRequest() {
		getContext().actorSelection(path).tell(new Identify(path), self());
	}

	@Override
	public void onReceive(Object msg) throws Exception {

		if (msg instanceof ActorIdentity) {
			master = ((ActorIdentity) msg).getRef();
			master.tell(WorkerMessage.IDLE, self());

		} else if (msg.equals(ReceiveTimeout.getInstance())) {
			sendIdentifyRequest();

		} else if (msg instanceof ControllerMessage) {
			if (msg == ControllerMessage.RUN_SIMULATION) {
				System.out.println("Asked to run simulation");
			}
		}

	}

}
