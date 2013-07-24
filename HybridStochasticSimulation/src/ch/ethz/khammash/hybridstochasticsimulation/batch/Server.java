package ch.ethz.khammash.hybridstochasticsimulation.batch;

import java.io.IOException;
import java.rmi.UnknownHostException;

import de.root1.simon.Registry;
import de.root1.simon.Simon;
import de.root1.simon.exceptions.NameBindingException;

public class Server {

	public static void main(String[] args) throws UnknownHostException, IOException, NameBindingException {

        // create the server's registry ...
        final Registry registry = Simon.createRegistry(22222);

        // create the serverobject
        ServerInterfaceImpl serverImpl = new ServerInterfaceImpl(
        		new ServerShutdownListener() {

					@Override
					public void shutdown() {
						registry.unbind("server");
						registry.stop();
					}

				}
        );

        // ... where we can bind the serverobject to
        registry.bind("server", serverImpl);

        System.out.println("Server up and running!");

        // some mechanism to shutdown the server should be placed here
        // this should include the following command:
        // registry.unbind("server");
        // registry.stop();
    }

}
