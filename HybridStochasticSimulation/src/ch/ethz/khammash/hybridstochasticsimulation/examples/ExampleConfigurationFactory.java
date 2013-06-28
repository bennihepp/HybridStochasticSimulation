package ch.ethz.khammash.hybridstochasticsimulation.examples;

import java.util.HashMap;
import java.util.Map;
import java.util.NoSuchElementException;
import java.util.Set;


public class ExampleConfigurationFactory {

	public interface ExampleConfigurationCreator {
		public ExampleConfiguration create();
	}

	private static ExampleConfigurationFactory instance = new ExampleConfigurationFactory();

	public static ExampleConfigurationFactory getInstance() {
		return instance;
	}

	private Map<String, ExampleConfigurationCreator> creatorMap;

	private ExampleConfigurationFactory() {
		creatorMap = new HashMap<String, ExampleConfigurationCreator>();
		registerExampleConfiguration("Trivial", new ExampleConfigurationCreator() {
			@Override
			public ExampleConfiguration create() {
				return new TrivialNetwork();
			}
		});
		registerExampleConfiguration("Simple Crystallization", new ExampleConfigurationCreator() {
			@Override
			public ExampleConfiguration create() {
				return new SimpleCrystallizationNetwork();
			}
		});
		registerExampleConfiguration("Birth Death Tunnel", new ExampleConfigurationCreator() {
			@Override
			public ExampleConfiguration create() {
				return new BirthDeathTunnelNetwork();
			}
		});
		registerExampleConfiguration("Haploinsufficiency", new ExampleConfigurationCreator() {
			@Override
			public ExampleConfiguration create() {
				return new HaploinsufficiencyNetwork();
			}
		});
		registerExampleConfiguration("Regulated Transcription", new ExampleConfigurationCreator() {
			@Override
			public ExampleConfiguration create() {
				return new RegulatedTranscriptionNetwork();
			}
		});
		registerExampleConfiguration("Bacterium Operator Site", new ExampleConfigurationCreator() {
			@Override
			public ExampleConfiguration create() {
				return new BacteriumOperatorSiteNetwork();
			}
		});
		registerExampleConfiguration("Lambda Phage Toggle Switch", new ExampleConfigurationCreator() {
			@Override
			public ExampleConfiguration create() {
				return new LambdaPhageToggleSwitchNetwork();
			}
		});
		registerExampleConfiguration("Repressed Bacterium Operon", new ExampleConfigurationCreator() {
			@Override
			public ExampleConfiguration create() {
				return new RepressedBacteriumOperonNetwork();
			}
		});
		registerExampleConfiguration("Conversion Cycle", new ExampleConfigurationCreator() {
			@Override
			public ExampleConfiguration create() {
				return new ConversionCycleNetwork();
			}
		});
		registerExampleConfiguration("Stochastic Focusing", new ExampleConfigurationCreator() {
			@Override
			public ExampleConfiguration create() {
				return new StochasticFocusingNetwork();
			}
		});
		registerExampleConfiguration("Heat Shock Response", new ExampleConfigurationCreator() {
			@Override
			public ExampleConfiguration create() {
				return new HeatShockResponseNetwork();
			}
		});
		registerExampleConfiguration("Vilar Oscillator", new ExampleConfigurationCreator() {
			@Override
			public ExampleConfiguration create() {
				return new VilarOscillator();
			}
		});
		registerExampleConfiguration("BacteriophageT7", new ExampleConfigurationCreator() {
			@Override
			public ExampleConfiguration create() {
				return new BacteriophageT7();
			}
		});
	}

	public void registerExampleConfiguration(String name, ExampleConfigurationCreator creator) {
		if (creatorMap.containsKey(name))
			throw new RuntimeException("Trying to register a duplicate ExampleConfiguration: \"" + name + "\"");
		creatorMap.put(name, creator);
	}

	public ExampleConfiguration createExampleConfiguration(String name) {
		if (!creatorMap.containsKey(name))
			throw new NoSuchElementException("No corresponding ExampleConfiguration has been registered: \"" + name + "\"");
		ExampleConfigurationCreator creator = creatorMap.get(name);
		return creator.create();
	}

	public Set<String> getAvailableExampleConfigurations() {
		return creatorMap.keySet();
	}

}
