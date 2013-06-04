package ch.ethz.khammash.hybridstochasticsimulation.examples;

import java.util.HashMap;
import java.util.Map;
import java.util.NoSuchElementException;
import java.util.Set;


public class ExampleNetworkFactory {

	public interface ExampleNetworkCreator {
		public ExampleNetwork create();
	}

	private static ExampleNetworkFactory instance = new ExampleNetworkFactory();

	public static ExampleNetworkFactory getInstance() {
		return instance;
	}

	private Map<String, ExampleNetworkCreator> creatorMap;

	private ExampleNetworkFactory() {
		creatorMap = new HashMap<String, ExampleNetworkCreator>();
		registerExampleNetwork("Simple Crystallization", new ExampleNetworkCreator() {
			@Override
			public ExampleNetwork create() {
				return new SimpleCrystallizationNetwork();
			}
		});
		registerExampleNetwork("Birth Death Tunnel", new ExampleNetworkCreator() {
			@Override
			public ExampleNetwork create() {
				return new BirthDeathTunnelNetwork();
			}
		});
		registerExampleNetwork("Haploinsufficiency", new ExampleNetworkCreator() {
			@Override
			public ExampleNetwork create() {
				return new HaploinsufficiencyNetwork();
			}
		});
		registerExampleNetwork("Regulated Transcription", new ExampleNetworkCreator() {
			@Override
			public ExampleNetwork create() {
				return new RegulatedTranscriptionNetwork();
			}
		});
		registerExampleNetwork("Bacterium Operator Site", new ExampleNetworkCreator() {
			@Override
			public ExampleNetwork create() {
				return new BacteriumOperatorSite();
			}
		});
		registerExampleNetwork("Lambda Phage Toggle Switch", new ExampleNetworkCreator() {
			@Override
			public ExampleNetwork create() {
				return new LambdaPhageToggleSwitchNetwork();
			}
		});
		registerExampleNetwork("Repressed Bacterium Operon", new ExampleNetworkCreator() {
			@Override
			public ExampleNetwork create() {
				return new RepressedBacteriumOperon();
			}
		});
	}

	public void registerExampleNetwork(String name, ExampleNetworkCreator creator) {
		if (creatorMap.containsKey(name))
			throw new RuntimeException("Trying to register a duplicate ExampleNetwork: \"" + name + "\"");
		creatorMap.put(name, creator);
	}

	public ExampleNetwork createExampleNetwork(String name) {
		if (!creatorMap.containsKey(name))
			throw new NoSuchElementException("No corresponding ExampleNetwork has been registered: \"" + name + "\"");
		ExampleNetworkCreator creator = creatorMap.get(name);
		return creator.create();
	}

	public Set<String> getAvailableExampleNetworks() {
		return creatorMap.keySet();
	}

}
