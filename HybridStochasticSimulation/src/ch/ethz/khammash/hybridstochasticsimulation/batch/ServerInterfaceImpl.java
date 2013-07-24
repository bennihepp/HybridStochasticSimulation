package ch.ethz.khammash.hybridstochasticsimulation.batch;

import de.root1.simon.annotation.SimonRemote;

@SimonRemote(value = {ServerInterface.class})
public class ServerInterfaceImpl implements ServerInterface {

	private ServerShutdownListener listener;

	ServerInterfaceImpl(ServerShutdownListener listener) {
		this.listener = listener;
	}

	@Override
	public void runSingleSimulation(ClientCallbackInterface callback) {
		callback.callback(null);
	}

	@Override
	public void shutdownServer() {
		listener.shutdown();
	}

}
