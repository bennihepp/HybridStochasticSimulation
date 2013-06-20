package ch.ethz.khammash.hybridstochasticsimulation;

import java.util.List;

import org.apache.commons.math3.util.FastMath;

import ch.ethz.khammash.hybridstochasticsimulation.gui.GUIEvent;
import ch.ethz.khammash.hybridstochasticsimulation.gui.GUIEvent.EventType;
import ch.ethz.khammash.hybridstochasticsimulation.gui.PlotWindow;
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
//							List<FinitePlotData> plotDataList = Main.trivialNetwork();
							// Maybe example
//							List<FinitePlotData> plotDataList = Main.conversionCycleNetwork();
//							List<FinitePlotData> plotDataList = Main.regulatedTranscriptionNetwork();
							// Example 2
//							List<FinitePlotData> plotDataList = Main.stochasticFocusingNetwork();
//							List<FinitePlotData> plotDataList = Main.simpleCrystallizationNetwork();
//							List<FinitePlotData> plotDataList = Main.birthDeathTunnelNetwork();
							// Example 1
//							List<FinitePlotData> plotDataList = Main.haploinsufficiencyNetwork();
//							List<FinitePlotData> plotDataList = Main.bacteriumOperatorSiteNetwork();
//							List<FinitePlotData> plotDataList = Main.lambdaPhageToggleSwitchNetwork();
//							List<FinitePlotData> plotDataList = Main.repressedBacteriumOperonNetwork();
							// Example 3
//							List<FinitePlotData> plotDataList = Main.heatShockResponseNetwork();
							List<FinitePlotData> plotDataList = Main.vilarOscillatorNetwork();
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
						}
					};
					t.start();
					break;
				case BENCHMARK:
					t = new Thread("Benchmark Thread") {
						@Override
						public void run() {
							window.getStatusBar().setText("Running benchmark");
							Main.vilarOscillatorNetwork();
							window.getStatusBar().setText("Finished benchmark", 2000, true);
						}
					};
					t.start();
					break;
				}
			}
		};
		window.getActionEventBus().register(actionHandler);
		window.pack();
		window.setVisible(true);
		window.getActionEventBus().post(new GUIEvent(EventType.SIMULATION));
	}

}
