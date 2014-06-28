package ch.ethz.bhepp.hybridstochasticsimulation.grid.akka;

import akka.actor.ActorSystem;
import akka.actor.Props;
import akka.kernel.Bootable;

import com.typesafe.config.Config;
import com.typesafe.config.ConfigFactory;
import com.typesafe.config.ConfigValueFactory;

public class AkkaWorkerSystem implements Bootable {

	private ActorSystem system;

	public AkkaWorkerSystem() {
		Config config = ConfigFactory.load("remote.conf").getConfig("worker");
		config = config.withValue("akka.remote.netty.tcp.port", ConfigValueFactory.fromAnyRef(2555));
		system = ActorSystem.create("AkkaGridWorker", config);
		final String path = "akka.tcp://AkkaGridController@127.0.0.1:2552/user/controller";
		system.actorOf(Props.create(AkkaWorker.class, path), "worker");
	}

	@Override
	public void startup() {
	}

	@Override
	public void shutdown() {
		system.shutdown();
	}

	public static void main(String[] args) {
		new AkkaWorkerSystem();
	}
}
