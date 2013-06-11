package ch.ethz.khammash.hybridstochasticsimulation.gui;

import static com.google.common.base.Preconditions.checkArgument;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Dimension;
import java.awt.GridLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.ComponentEvent;
import java.awt.event.ComponentListener;
import java.awt.event.InputEvent;
import java.io.File;
import java.io.IOException;
import java.util.Collection;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;

import javax.swing.JFileChooser;
import javax.swing.JFrame;
import javax.swing.JMenu;
import javax.swing.JMenuBar;
import javax.swing.JMenuItem;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.KeyStroke;
import javax.swing.filechooser.FileNameExtensionFilter;

import matlabcontrol.MatlabConnectionException;
import matlabcontrol.MatlabInvocationException;

import org.jfree.chart.ChartPanel;
import org.jfree.ui.ApplicationFrame;

import ch.ethz.khammash.hybridstochasticsimulation.matlab.MatlabDataExporter;
import ch.ethz.khammash.hybridstochasticsimulation.matlab.MatlabPlotter;
import ch.ethz.khammash.hybridstochasticsimulation.matlab.MatlabSession;
import ch.ethz.khammash.hybridstochasticsimulation.trajectories.PlotData;
import ch.ethz.khammash.hybridstochasticsimulation.trajectories.TrajectoryDistributionPlotData;
import ch.ethz.khammash.hybridstochasticsimulation.trajectories.TrajectoryPlotData;

import com.google.common.base.Optional;
import com.google.common.eventbus.EventBus;
import com.jmatio.types.MLArray;


public class PlotWindow extends ApplicationFrame implements ComponentListener {
	private static final long serialVersionUID = 6566839328303930162L;

	public static enum EventType {
		SIMULATION, BENCHMARK,
	}

	public static class Event {

		private EventType type;
		private ActionEvent event;

		public Event(EventType type) {
			this(type, null);
		}

		public Event(EventType type, ActionEvent event) {
			this.type = type;
			this.event = event;
		}

		public EventType getDescription() {
			return type;
		}

		public ActionEvent getEvent() {
			return event;
		}

	}

	private EventBus actionEventBus;

	private JScrollPane scrollPane;
	private JPanel mainPanel;
	private List<PlotData> plotDataList;
	private Map<PlotData, Optional<ChartPanel>> plotDataMap;
	private int rows;
	private int cols;

	private int savePlotWidth = 1600;
	private int savePlotHeight = 800;
	private int saveHorizontalSpacing = 10;
	private int saveVerticalSpacing = 10;

	public PlotWindow(String title) {
		this(title, 1, 1);
	}

	public PlotWindow(String title, int rows, int cols) {
		super(title);
		actionEventBus = new EventBus("Action Event Bus");
		checkArgument(rows > 0);
		checkArgument(cols > 0);
		plotDataList = new LinkedList<PlotData>();
		plotDataMap = new LinkedHashMap<PlotData, Optional<ChartPanel>>();
		setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		setJMenuBar(createMenuBar());
		setPreferredSize(new Dimension(1400, 700));
		mainPanel = new PlotPanel();
		scrollPane = new JScrollPane(mainPanel,
				JScrollPane.VERTICAL_SCROLLBAR_AS_NEEDED,
				JScrollPane.HORIZONTAL_SCROLLBAR_AS_NEEDED);
		scrollPane.getViewport().setBackground(Color.WHITE);
		setLayout(new BorderLayout());
		add(scrollPane, BorderLayout.CENTER);
		updateLayout(rows, cols);
	}

	final private void updateLayout(int rows, int cols) {
		clearPlotData();
		this.rows = rows;
		this.cols = cols;
		mainPanel.setLayout(new GridLayout(rows, cols));
	}

	public void clearPlotData() {
		Iterator<PlotData> it = plotDataList.iterator();
		while (it.hasNext()) {
			PlotData plotData = it.next();
			it.remove();
			Optional<ChartPanel> optionalPanel = plotDataMap.remove(plotData);
			if (optionalPanel.isPresent())
				mainPanel.remove(optionalPanel.get());
		}
	}

	public void addPlotData(Collection<PlotData> plotDataCollection) {
		Iterator<PlotData> it = plotDataCollection.iterator();
		while (it.hasNext())
			addPlotData(it.next());
	}

	public void setPlotData(Collection<PlotData> plotDataCollection, int rows, int cols) {
		updateLayout(rows, cols);
		addPlotData(plotDataCollection);
	}

	public void addPlotData(PlotData plotData) {
		Optional<ChartPanel> optionalPanel = Optional.<ChartPanel>absent();
		if (plotData instanceof TrajectoryDistributionPlotData) {
			TrajectoryDistributionPlotChartPanel panel = new TrajectoryDistributionPlotChartPanel();
			panel.addDistributionPlotData((TrajectoryDistributionPlotData)plotData);
			mainPanel.add(panel);
			optionalPanel = Optional.<ChartPanel>of(panel);
		} else if (plotData instanceof TrajectoryPlotData) {
			TrajectoryPlotChartPanel panel = new TrajectoryPlotChartPanel();
			panel.addPlotData((TrajectoryPlotData)plotData);
			mainPanel.add(panel);
			optionalPanel = Optional.<ChartPanel>of(panel);
		} else
			throw new UnsupportedOperationException("Unsupported type of PlotData");
		plotDataMap.put(plotData, optionalPanel);
		plotDataList.add(plotData);
	}

	protected JMenuBar createMenuBar() {
		// Export to Matlab MAT file entry
		JMenuItem exportToMatlabFileEntry = new JMenuItem("Export to Matlab MAT file");
		exportToMatlabFileEntry.addActionListener(
				new ActionListener() {
					@Override
					public void actionPerformed(ActionEvent e) {
						JFileChooser fc = new CustomJFileChooser();
						fc.setFileFilter(new FileNameExtensionFilter("MAT files", "mat"));
						int returnValue = fc.showSaveDialog(PlotWindow.this);
						if (returnValue == JFileChooser.APPROVE_OPTION)
							exportToMatlabFile(fc.getSelectedFile());
					}
				}
			);
		// Plot in Matlab entry
		JMenuItem plotInMatlabEntry = new JMenuItem("Plot in Matlab");
		plotInMatlabEntry.setAccelerator(KeyStroke.getKeyStroke('M', InputEvent.CTRL_MASK));
		plotInMatlabEntry.addActionListener(
				new ActionListener() {
					@Override
					public void actionPerformed(ActionEvent e) {
						plotInMatlab();
					}
				}
			);
		// Quit entry
		JMenuItem quitEntry = new JMenuItem("Quit");
		quitEntry.setAccelerator(KeyStroke.getKeyStroke('Q', InputEvent.CTRL_MASK));
		quitEntry.addActionListener(
			new ActionListener() {
				@Override
				public void actionPerformed(ActionEvent e) {
					PlotWindow.this.dispose();
				}
			}
		);
		// Run entry
		JMenuItem runEntry = new JMenuItem("Run");
		runEntry.setAccelerator(KeyStroke.getKeyStroke('R', InputEvent.CTRL_MASK));
		runEntry.addActionListener(
			new ActionListener() {
				@Override
				public void actionPerformed(ActionEvent e) {
					getActionEventBus().post(new Event(EventType.SIMULATION, e));
				}
			}
		);
		// Benchmark entry
		JMenuItem bechmarkEntry = new JMenuItem("Benchmark");
		bechmarkEntry.setAccelerator(KeyStroke.getKeyStroke('B', InputEvent.CTRL_MASK));
		bechmarkEntry.addActionListener(
			new ActionListener() {
				@Override
				public void actionPerformed(ActionEvent e) {
					getActionEventBus().post(new Event(EventType.BENCHMARK, e));
				}
			}
		);
		// Add entries to menubar
		JMenuBar menubar = new JMenuBar();
		JMenu fileMenu = new JMenu("File");
		fileMenu.add(exportToMatlabFileEntry);
		fileMenu.add(plotInMatlabEntry);
		fileMenu.addSeparator();
		fileMenu.add(quitEntry);
		JMenu simulationMenu = new JMenu("Simulation");
		simulationMenu.add(runEntry);
		simulationMenu.add(bechmarkEntry);
		menubar.add(fileMenu);
		menubar.add(simulationMenu);
		return menubar;
	}

	public void exportToMatlabFile(File file) {
		MatlabDataExporter mde = new MatlabDataExporter();
		List<MLArray> matlabData = mde.buildMatlabData(plotDataList, rows, cols);
		try {
			mde.writeMatlabDataToFile(file, matlabData);
		} catch (IOException e1) {
			String errorMsg = "Failed to export Matlab file:\n" + e1.toString();
			e1.printStackTrace();
			JOptionPane.showMessageDialog(PlotWindow.this, errorMsg, "Error", JOptionPane.ERROR_MESSAGE);
		}
	}

	public void plotInMatlab() {
		MatlabSession session = new MatlabSession();
		try {
			session.start();
			MatlabPlotter mp = new MatlabPlotter(session);
			mp.plot(plotDataList, rows, cols);
		} catch (MatlabConnectionException e1) {
			String errorMsg = "Failed to plot in Matlab:\n" + e1.toString();
			e1.printStackTrace();
			JOptionPane.showMessageDialog(PlotWindow.this, errorMsg, "Error", JOptionPane.ERROR_MESSAGE);
		} catch (MatlabInvocationException e1) {
			String errorMsg = "Failed to plot in Matlab:\n" + e1.toString();
			e1.printStackTrace();
			JOptionPane.showMessageDialog(PlotWindow.this, errorMsg, "Error", JOptionPane.ERROR_MESSAGE);
		}
	}

	public int getSavePlotWidth() {
		return savePlotWidth;
	}

	public void setSavePlotWidth(int savePlotWidth) {
		this.savePlotWidth = savePlotWidth;
	}

	public int getSavePlotHeight() {
		return savePlotHeight;
	}

	public void setSavePlotHeight(int savePlotHeight) {
		this.savePlotHeight = savePlotHeight;
	}

	public int getSaveHorizontalSpacing() {
		return saveHorizontalSpacing;
	}

	public void setSaveHorizontalSpacing(int saveHorizontalSpacing) {
		this.saveHorizontalSpacing = saveHorizontalSpacing;
	}

	public int getSaveVerticalSpacing() {
		return saveVerticalSpacing;
	}

	public void setSaveVerticalSpacing(int saveVerticalSpacing) {
		this.saveVerticalSpacing = saveVerticalSpacing;
	}

	public EventBus getActionEventBus() {
		return actionEventBus;
	}

	@Override
	public void componentResized(ComponentEvent e) {
//		System.out.println("resized");
//		Dimension scrollPaneSize = scrollPane.getViewport().getExtentSize();
//		Dimension mainPanelSize = mainPanel.getSize();
//		mainPanelSize.setSize(scrollPaneSize.getWidth(), mainPanelSize.getHeight());
//		mainPanel.setPreferredSize(scrollPane.getViewport().getExtentSize());
	}

	@Override
	public void componentMoved(ComponentEvent e) {
	}

	@Override
	public void componentShown(ComponentEvent e) {
	}

	@Override
	public void componentHidden(ComponentEvent e) {
	}

}
