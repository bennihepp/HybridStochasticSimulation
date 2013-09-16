package ch.ethz.khammash.hybridstochasticsimulation.grid.akka;

import akka.actor.ActorSystem;
import akka.actor.Props;
import akka.kernel.Bootable;

import com.typesafe.config.ConfigFactory;

public class AkkaControllerSystem implements Bootable {

	private ActorSystem system;

	public AkkaControllerSystem() {
		system = ActorSystem.create("AkkaGridController", ConfigFactory.load().getConfig("controller"));
		system.actorOf(Props.create(AkkaController.class), "controller");
	}

	@Override
	public void startup() {
	}

	@Override
	public void shutdown() {
		system.shutdown();
	}

	public static void main(String[] args) {
		new AkkaControllerSystem();
	}

}
