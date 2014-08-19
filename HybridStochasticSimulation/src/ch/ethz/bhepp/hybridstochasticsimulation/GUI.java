package ch.ethz.bhepp.hybridstochasticsimulation;

import java.util.List;

import org.apache.commons.math3.util.FastMath;

import ch.ethz.bhepp.hybridstochasticsimulation.gui.GUIEvent;
import ch.ethz.bhepp.hybridstochasticsimulation.gui.PlotWindow;
import ch.ethz.bhepp.hybridstochasticsimulation.gui.GUIEvent.EventType;
import ch.ethz.bhepp.hybridstochasticsimulation.trajectories.FinitePlotData;

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
							try {
//								List<FinitePlotData> plotDataList = Examples.fastIsomerization();
//								List<FinitePlotData> plotDataList = Examples.fastDimerization();
								List<FinitePlotData> plotDataList = Examples.bacteriophageT7Network();
//								List<FinitePlotData> plotDataList = Examples.repressilator();
//								List<FinitePlotData> plotDataList = Examples.toggleSwitch();
//								List<FinitePlotData> plotDataList = Examples.haploinsufficiencyNetwork();								int rows;

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
							} catch (InterruptedException e) {
								window.getStatusBar().setText("Simulation was interrupted", 2000);
							}
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
