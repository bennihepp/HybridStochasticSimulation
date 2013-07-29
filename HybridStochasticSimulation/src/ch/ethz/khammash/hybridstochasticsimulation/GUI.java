package ch.ethz.khammash.hybridstochasticsimulation;

import java.io.IOException;
import java.util.List;

import javax.xml.parsers.ParserConfigurationException;

import org.apache.commons.math3.util.FastMath;
import org.xml.sax.SAXException;

import ch.ethz.khammash.hybridstochasticsimulation.gui.GUIEvent;
import ch.ethz.khammash.hybridstochasticsimulation.gui.GUIEvent.EventType;
import ch.ethz.khammash.hybridstochasticsimulation.gui.PlotWindow;
import ch.ethz.khammash.hybridstochasticsimulation.io.StochKitNetworkReader.FileFormatException;
import ch.ethz.khammash.hybridstochasticsimulation.trajectories.FinitePlotData;

import com.google.common.eventbus.Subscribe;

public class GUI {

	public static void main(String[] args) {
		final PlotWindow window = new PlotWindow("Hybrid stochastic simulations");
		Object actionHandler = new Object() {
			@Subscribe
			public void handleActionEvent(GUIEvent e) {
				switch (e.getDescription()) {
				case SIMULATION:
					Thread t = new Thread("Simulation Thread") {
						@Override
						public void run() {
							window.getStatusBar().setText("Running simulation");
//							List<FinitePlotData> plotDataList = Examples.trivialNetwork();
							// Maybe example
//							List<FinitePlotData> plotDataList = Examples.regulatedTranscriptionNetwork();
							// Example 2
//							List<FinitePlotData> plotDataList = Examples.birthDeathTunnelNetwork();
							// Example 1
//							List<FinitePlotData> plotDataList = Examples.bacteriumOperatorSiteNetwork();
//							List<FinitePlotData> plotDataList = Examples.lambdaPhageToggleSwitchNetwork();
//							List<FinitePlotData> plotDataList = Examples.repressedBacteriumOperonNetwork();
							// Example 3
//							List<FinitePlotData> plotDataList = Examples.heatShockResponseNetwork();
//							List<FinitePlotData> plotDataList = Examples.vilarOscillatorNetwork();

//							List<FinitePlotData> plotDataList = Examples.bacteriophageT7Network();
//							List<FinitePlotData> plotDataList = Examples.fastIsomerization();
//							List<FinitePlotData> plotDataList = Examples.fastDimerization();
//							List<FinitePlotData> plotDataList = Examples.haploinsufficiencyNetwork();
//							List<FinitePlotData> plotDataList = Examples.simpleCrystallizationNetwork();
//							List<FinitePlotData> plotDataList = Examples.stochasticFocusingNetwork();
//							List<FinitePlotData> plotDataList = Examples.conversionCycleNetwork();

//							List<FinitePlotData> plotDataList = Examples.repressilator();
//							List<FinitePlotData> plotDataList = Examples.toggleSwitch();
							try {
//								List<FinitePlotData> plotDataList = Examples.heatShockMassAction();
//								List<FinitePlotData> plotDataList = Examples.enzymeKinetics1();
								List<FinitePlotData> plotDataList = Examples.complexExample1();
								int rows;
								int cols;
								if (plotDataList.size() >= 3) {
									rows = (int) FastMath.ceil(plotDataList.size() / 3.0);
									cols = 3;
								} else {
									rows = 1;
									cols = (int) FastMath.ceil(plotDataList.size() / 3.0);
								}
								window.clearPlotData();
								window.setPlotData(plotDataList, rows, cols);
								window.validate();
								window.getStatusBar().setText("Finished simulation", 2000);
							} catch (ParserConfigurationException
									| SAXException | IOException | FileFormatException e) {
								e.printStackTrace();
								window.getStatusBar().setText("Simulation failed!");
							}
						}
					};
					t.start();
					break;
				case BENCHMARK:
					t = new Thread("Benchmark Thread") {
						@Override
						public void run() {
							window.getStatusBar().setText("Running benchmark");
							Examples.vilarOscillatorNetwork();
							window.getStatusBar().setText("Finished benchmark", 2000, true);
						}
					};
					t.start();
					break;
				}
			}
		};

//		ExampleConfiguration nss = ExampleConfigurationFactory.getInstance().createExampleConfiguration("Vilar Oscillator");
//
//    	ReactionNetworkGraph graph = new ReactionNetworkGraph(nss.net);
//    	DependencyGraph depGraph = new DependencyGraph(graph);
//    	for (SpeciesVertex v : depGraph.vertexSet()) {
//    		System.out.println("Reachable from " + v);
//    		for (DependencyEdge e : depGraph.outgoingEdgesOf(v)) {
//    			System.out.println("  " + e.getTarget());
//    		}
//    	}

//        // create a JGraphT graph
//    	ReactionNetworkGraph graph = new ReactionNetworkGraph(nss.net);
//    	HashMap<Integer,SpeciesVertex> vertexMap = new HashMap<Integer, SpeciesVertex>();
//    	for (SpeciesVertex v : graph.vertexSet()) {
//    		v.setName(nss.speciesNames[v.getSpecies()]);
//    		vertexMap.put(v.getSpecies(), v);
//    	}
//
//		GraphWindow gwindow = new GraphWindow(graph);
//		gwindow.pack();
////		gwindow.autoPosition();
//		double width = gwindow.getWidth();
//		double height = gwindow.getHeight();
//		gwindow.positionVertexAt(vertexMap.get(0), width * 4 / 8., height * 4 / 8.);
//		gwindow.positionVertexAt(vertexMap.get(1), width * 3 / 4., height * 7 / 8.);
//		gwindow.positionVertexAt(vertexMap.get(2), width * 5 / 6., height * 2 / 6.);
//		gwindow.positionVertexAt(vertexMap.get(3), width * 3 / 6., height * 7 / 8.);
//		gwindow.positionVertexAt(vertexMap.get(4), width * 5 / 6., height * 4 / 6.);
//		gwindow.positionVertexAt(vertexMap.get(5), width * 1 / 6., height * 2 / 6.);
//		gwindow.positionVertexAt(vertexMap.get(6), width * 2 / 6., height * 1 / 6.);
//		gwindow.positionVertexAt(vertexMap.get(7), width * 3 / 6., height * 1 / 6.);
//		gwindow.positionVertexAt(vertexMap.get(8), width * 1 / 6., height * 5 / 8.);
////		gwindow.revalidate();
//		gwindow.setVisible(true);

		window.getActionEventBus().register(actionHandler);
		window.pack();
		window.setVisible(true);
		window.getActionEventBus().post(new GUIEvent(EventType.SIMULATION));

	}

}
