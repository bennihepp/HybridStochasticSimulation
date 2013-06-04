package ch.ethz.khammash.hybridstochasticsimulation.plotting;

import static com.google.common.base.Preconditions.checkArgument;

import java.awt.Component;
import java.awt.GridLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.InputEvent;
import java.io.IOException;
import java.util.LinkedList;
import java.util.List;

import javax.swing.JFileChooser;
import javax.swing.JFrame;
import javax.swing.JMenu;
import javax.swing.JMenuBar;
import javax.swing.JMenuItem;
import javax.swing.JOptionPane;
import javax.swing.KeyStroke;
import javax.swing.filechooser.FileNameExtensionFilter;

import matlabcontrol.MatlabConnectionException;
import matlabcontrol.MatlabInvocationException;

import org.jfree.ui.ApplicationFrame;

import com.jmatio.types.MLArray;


public class GridWindow extends ApplicationFrame {
	private static final long serialVersionUID = 6566839328303930162L;

	List<PlotData> plotDataList;
	int rows;
	int cols;

	private int savePlotWidth = 1600;
	private int savePlotHeight = 800;
	private int saveHorizontalSpacing = 10;
	private int saveVerticalSpacing = 10;

	public GridWindow(String title, int rows, int cols) {
		super(title);

		checkArgument(rows > 0);
		checkArgument(cols > 0);
		plotDataList = new LinkedList<PlotData>();
		this.rows = rows;
		this.cols = cols;
		setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		setLayout(new GridLayout(rows, cols));
		setJMenuBar(createMenuBar());
	}

	public void addPlotData(PlotData plotData) {
		if (plotData instanceof TrajectoryDistributionPlotData) {
			TrajectoryDistributionPlotChartPanel panel = new TrajectoryDistributionPlotChartPanel();
			panel.addDistributionPlotData((TrajectoryDistributionPlotData)plotData);
			add(panel);
		} else if (plotData instanceof TrajectoryPlotData) {
			TrajectoryPlotChartPanel panel = new TrajectoryPlotChartPanel();
			panel.addPlotData((TrajectoryPlotData)plotData);
			add(panel);
		} else
			throw new UnsupportedOperationException("Unsupported type of PlotData");
		plotDataList.add(plotData);
	}

	@Override
	public Component add(Component comp) {
		Component q = super.add(comp);
		pack();
		return q;
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
						int returnValue = fc.showSaveDialog(GridWindow.this);
						if (returnValue == JFileChooser.APPROVE_OPTION) {
							MatlabDataExporter mde = new MatlabDataExporter();
							List<MLArray> matlabData = mde.buildMatlabData(plotDataList, rows, cols);
							try {
								mde.writeMatlabDataToFile(fc.getSelectedFile(), matlabData);
							} catch (IOException e1) {
								String errorMsg = "Failed to export Matlab file:\n" + e1.toString();
								JOptionPane.showMessageDialog(GridWindow.this, errorMsg, "Error", JOptionPane.ERROR);
							}
						}
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
						MatlabSession session = new MatlabSession();
						try {
							session.start();
							MatlabPlotter mp = new MatlabPlotter(session);
							mp.plot(plotDataList, rows, cols);
						} catch (MatlabConnectionException e1) {
							String errorMsg = "Failed to plot in Matlab:\n" + e1.toString();
							JOptionPane.showMessageDialog(GridWindow.this, errorMsg, "Error", ERROR);
						} catch (MatlabInvocationException e1) {
							String errorMsg = "Failed to plot in Matlab:\n" + e1.toString();
							JOptionPane.showMessageDialog(GridWindow.this, errorMsg, "Error", ERROR);
						}
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
					GridWindow.this.dispose();
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
		menubar.add(fileMenu);
		return menubar;
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
}
