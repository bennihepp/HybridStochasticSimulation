package ch.ethz.bhepp.hybridstochasticsimulation.examples;

import java.util.HashMap;
import java.util.Map;
import java.util.NoSuchElementException;
import java.util.Set;


public class ExampleConfigurationFactory {

	public interface ExampleConfigurationCreator {
		public SimulationConfiguration create();
	}

	private static volatile ExampleConfigurationFactory instance = null;

	public static ExampleConfigurationFactory getInstance() {
		if (instance == null) {
			synchronized (ExampleConfigurationFactory.class) {
				if (instance == null) {
					instance = new ExampleConfigurationFactory();
				}
			}
		}
		return instance;
	}

	private Map<String, ExampleConfigurationCreator> creatorMap;

	private ExampleConfigurationFactory() {
		creatorMap = new HashMap<String, ExampleConfigurationCreator>();
		registerExampleConfiguration("Trivial", new ExampleConfigurationCreator() {
			@Override
			public SimulationConfiguration create() {
				return new TrivialNetwork();
			}
		});
		registerExampleConfiguration("Simple Crystallization", new ExampleConfigurationCreator() {
			@Override
			public SimulationConfiguration create() {
				return new SimpleCrystallizationNetwork();
			}
		});
		registerExampleConfiguration("Birth Death Tunnel", new ExampleConfigurationCreator() {
			@Override
			public SimulationConfiguration create() {
				return new BirthDeathTunnelNetwork();
			}
		});
		registerExampleConfiguration("Haploinsufficiency", new ExampleConfigurationCreator() {
			@Override
			public SimulationConfiguration create() {
				return new HaploinsufficiencyNetwork();
			}
		});
		registerExampleConfiguration("Regulated Transcription", new ExampleConfigurationCreator() {
			@Override
			public SimulationConfiguration create() {
				return new RegulatedTranscriptionNetwork();
			}
		});
		registerExampleConfiguration("Bacterium Operator Site", new ExampleConfigurationCreator() {
			@Override
			public SimulationConfiguration create() {
				return new BacteriumOperatorSiteNetwork();
			}
		});
		registerExampleConfiguration("Lambda Phage Toggle Switch", new ExampleConfigurationCreator() {
			@Override
			public SimulationConfiguration create() {
				return new LambdaPhageToggleSwitchNetwork();
			}
		});
		registerExampleConfiguration("Repressed Bacterium Operon", new ExampleConfigurationCreator() {
			@Override
			public SimulationConfiguration create() {
				return new RepressedBacteriumOperonNetwork();
			}
		});
		registerExampleConfiguration("Conversion Cycle", new ExampleConfigurationCreator() {
			@Override
			public SimulationConfiguration create() {
				return new ConversionCycleNetwork();
			}
		});
		registerExampleConfiguration("Stochastic Focusing", new ExampleConfigurationCreator() {
			@Override
			public SimulationConfiguration create() {
				return new StochasticFocusingNetwork();
			}
		});
		registerExampleConfiguration("Heat Shock Response", new ExampleConfigurationCreator() {
			@Override
			public SimulationConfiguration create() {
				return new HeatShockResponseNetwork();
			}
		});
		registerExampleConfiguration("Vilar Oscillator", new ExampleConfigurationCreator() {
			@Override
			public SimulationConfiguration create() {
				return new VilarOscillator();
			}
		});
		registerExampleConfiguration("BacteriophageT7", new ExampleConfigurationCreator() {
			@Override
			public SimulationConfiguration create() {
				return new BacteriophageT7();
			}
		});
		registerExampleConfiguration("FastIsomerization", new ExampleConfigurationCreator() {
			@Override
			public SimulationConfiguration create() {
				return new FastIsomerization();
			}
		});
		registerExampleConfiguration("FastDimerization", new ExampleConfigurationCreator() {
			@Override
			public SimulationConfiguration create() {
				return new FastDimerization();
			}
		});
		registerExampleConfiguration("ExtendedFastDimerization", new ExampleConfigurationCreator() {
			@Override
			public SimulationConfiguration create() {
				return new ExtendedFastDimerization();
			}
		});
		registerExampleConfiguration("ToggleSwitch", new ExampleConfigurationCreator() {
			@Override
			public SimulationConfiguration create() {
				return new ToggleSwitch();
			}
		});
		registerExampleConfiguration("Repressilator", new ExampleConfigurationCreator() {
			@Override
			public SimulationConfiguration create() {
				return new Repressilator();
			}
		});
		registerExampleConfiguration("ModifiedRepressilator", new ExampleConfigurationCreator() {
			@Override
			public SimulationConfiguration create() {
				return new ModifiedRepressilator();
			}
		});
	}

	public void registerExampleConfiguration(String name, ExampleConfigurationCreator creator) {
		if (creatorMap.containsKey(name))
			throw new RuntimeException("Trying to register a duplicate ExampleConfiguration: \"" + name + "\"");
		creatorMap.put(name, creator);
	}

	public SimulationConfiguration createExampleConfiguration(String name) {
		if (!creatorMap.containsKey(name))
			throw new NoSuchElementException("No corresponding ExampleConfiguration has been registered: \"" + name + "\"");
		ExampleConfigurationCreator creator = creatorMap.get(name);
		return creator.create();
	}

	public Set<String> getAvailableExampleConfigurations() {
		return creatorMap.keySet();
	}

}
