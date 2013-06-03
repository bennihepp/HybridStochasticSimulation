package ch.ethz.khammash.hybridstochasticsimulation.plotting;

import static com.google.common.base.Preconditions.checkArgument;

import java.awt.Component;
import java.awt.GridLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.InputEvent;
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

//	public int print(Graphics g, PageFormat pf, int page) throws PrinterException {
//	    if (page > 0) {
//	        return NO_SUCH_PAGE;
//	    }
//	    double x0 = pf.getImageableX();
//	    double y0 = pf.getImageableY();
//	    double width = pf.getImageableWidth();
//	    double height = pf.getImageableHeight();
//	    double plotWidth = width / rows;
//	    double plotHeight = height / rows;
//	    Graphics2D g2d = (Graphics2D)g;
//		int row = 0;
//		int col = 0;
//		for (ChartPanel panel : GridWindow.this.chartPanels) {
//			JFreeChart chart = panel.getChart();
//			double x = x0 + plotWidth * col;
//			double y = y0 + plotHeight * row;
//			java.awt.Rectangle area = new java.awt.Rectangle((int)x, (int)y, (int)plotWidth, (int)plotHeight);
//			chart.draw(g2d, area);
//			col++;
//			if (col >= cols) {
//				col = 0;
//				row++;
//			}
//		}
//	    return PAGE_EXISTS;
//	}

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
							mde.writeMatlabDataToFile(fc.getSelectedFile(), matlabData);
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
//		// Save as PDF with R entry
//		JMenuItem saveAsPdfWithREntry = new JMenuItem("Save as PDF with R");
//		saveAsPdfWithREntry.addActionListener(
//				new ActionListener() {
//					@Override
//					public void actionPerformed(ActionEvent e) {
//						JFileChooser fc = new CustomJFileChooser();
//						fc.setFileFilter(new FileNameExtensionFilter("PDF documents", "pdf"));
//						int returnValue = fc.showSaveDialog(GridWindow.this);
//						if (returnValue == JFileChooser.APPROVE_OPTION) {
//							// Just making sure we have the right version of everything
//							if (!Rengine.versionCheck())
//							    throw new UnsupportedClassVersionError("JRI version mismatch");
//							String[] args = { "--no-save", "--vanilla" };
//							Rengine re=new Rengine(args, false, new RConsole(GridWindow.this));
//					        if (!re.waitForR())
//							    throw new RuntimeException("Cannot load R");
////					        re.eval(".setenv <- if (exists(\"Sys.setenv\")) Sys.setenv else Sys.putenv");
////					        re.eval("library(JavaGD)");
////					        re.eval("JavaGD()");
//							File outputFile = fc.getSelectedFile();
//					        re.eval("pdf(\"" + outputFile.getAbsolutePath() + "\", width=11.7, height=8.3, pointsize=" + (12.0 / 600 * 72) + ")");
//					        re.eval("par(mfrow=c(" + rows + "," + cols + "))");
//							for (ChartPanel panel : chartPanels) {
//						        re.eval("par(col=\"black\")");
//								JFreeChart chart = panel.getChart();
//								XYPlot plot = (XYPlot)chart.getPlot();
//								XYSeriesCollection seriesCollection = (XYSeriesCollection)plot.getDataset();
//						        re.eval("colors <- rainbow(" + seriesCollection.getSeriesCount() + ")");
//								for (int j=0; j < seriesCollection.getSeriesCount(); j++) {
//									XYSeries series = (XYSeries)seriesCollection.getSeries(j);
//									double[] tSeries = new double[series.getItemCount()];
//									double[] xSeries = new double[series.getItemCount()];
//									for (int k=0; k < series.getItemCount(); k++) {
//										XYDataItem dataItem = series.getDataItem(k);
//										tSeries[k] = dataItem.getXValue();
//										xSeries[k] = dataItem.getYValue();
//									}
//									re.assign("tSeries", tSeries);
//									re.assign("xSeries", xSeries);
//									if (j == 0)
//								        re.eval("plot(tSeries, xSeries, type=\"n\", xlab=\"" + plot.getDomainAxis().getLabel() + "\", ylab=\"" + plot.getRangeAxis().getLabel() + "\")");
//									re.eval("lines(tSeries, xSeries, type=\"l\", col=colors[" + (j+1) + "])");
//								}
//								re.eval("title(\"" + chart.getTitle().getText() + "\")");
//								String[] plotLabels = new String[plot.getLegendItems().getItemCount()];
//								for (int j=0; j < plot.getLegendItems().getItemCount(); j++)
//									plotLabels[j] = plot.getLegendItems().get(j).getLabel();
//								re.assign("plotLabels",  plotLabels);
//								re.eval("legend(\"topright\", NULL, plotLabels, colors)");
//							}
//							re.eval("dev.off()");
//						    re.end();
//						}
//					}
//				}
//			);
//		// Test entry
//		JMenuItem testEntry = new JMenuItem("Test");
//		testEntry.addActionListener(
//				new ActionListener() {
//					@Override
//					public void actionPerformed(ActionEvent e) {
//						String fileName = "C:\\Users\\Benni\\Documents\\test3.pdf";
//						int width = 3200;
//						int height = 2400;
//						JFreeChart chart = chartPanels.get(0).getChart();
//					    if (chart != null) {
//					        BufferedOutputStream out = null;
//					        try {
//					            out = new BufferedOutputStream(new FileOutputStream(fileName)); 
//
//					            //convert chart to PDF with iText:
//					            Rectangle pagesize = new Rectangle(width, height); 
//					            Document document = new Document(pagesize, 50, 50, 50, 50); 
//					            try { 
//					                PdfWriter writer = PdfWriter.getInstance(document, out); 
//					                document.addAuthor("JFreeChart"); 
//					                document.open(); 
//
//					                DefaultFontMapper fontMapper = new DefaultFontMapper();
//					                fontMapper.insertDirectory("resources/");
//					                Font font = UIManager.getDefaults().getFont("TabbedPane.font");
//					                GraphicsEnvironment.getLocalGraphicsEnvironment().registerFont(font);
//
//									StandardChartTheme theme = new StandardChartTheme("CustomChartFont", true);
//								    theme.setExtraLargeFont(font.deriveFont(Font.BOLD, 20*5));
//								    theme.setLargeFont(font.deriveFont(Font.BOLD, 16*5));
//								    theme.setRegularFont(font.deriveFont(Font.PLAIN, 14*5));
//								    theme.setSmallFont(font.deriveFont(Font.PLAIN, 12*5));
//								    ChartFactory.setChartTheme(theme);
//
//					                PdfContentByte cb = writer.getDirectContent(); 
////					                PdfTemplate tp = cb.createTemplate(width, height); 
//					                Graphics2D g2 = new PdfGraphics2D(cb, width, height, fontMapper);
////					                Graphics2D g2 = tp.createGraphics(width, height, new DefaultFontMapper()); 
//
//					                Rectangle2D r2D = new Rectangle2D.Double(0, 0, width, height); 
//					                chart.draw(g2, r2D, null); 
//					                g2.dispose(); 
////					                cb.addTemplate(tp, 0, 0);
//					            } catch (DocumentException e1) {
//					            	e1.printStackTrace();
//					            } finally {
//					                document.close(); 
//					            }
//				            } catch (IOException e1) {
//				            	e1.printStackTrace();
//					        } finally {
//					            if (out != null) 
//					            	try {
//					            		out.close();
//					            	} catch (IOException e1) {
//					            	}
//					        }
//					    }//else: input values not availabel
//					}
//				}
//			);
//		// Print entry
//		JMenuItem printEntry = new JMenuItem("Print");
//		printEntry.addActionListener(
//				new ActionListener() {
//					@Override
//					public void actionPerformed(ActionEvent e) {
//						final int DPI = 300;
////						PrinterResolution pr = new PrinterResolution(DPI, DPI, PrinterResolution.DPI);
////					    PrintRequestAttributeSet aset = new HashPrintRequestAttributeSet();
////					    aset.add(MediaSizeName.ISO_A4);
////					    aset.add(pr);
////					    aset.add(OrientationRequested.PORTRAIT);
////
////					    Toolkit tk = Toolkit.getDefaultToolkit();
////					    PrintJob pj = tk.getPrintJob(GridWindow.this, "Print plots", null);
////					    Graphics g = pj.getGraphics();
//
//						final double a4WidthInch = 8.26771654; // Equals 210mm
//		                final double a4HeightInch = 11.6929134; // Equals 297mm
//		                final int orientation = PageFormat.LANDSCAPE;
//					    Paper paper = new Paper();
//					    paper.setSize(72.0 * a4WidthInch * DPI / 72.0, 72.0 * a4HeightInch * DPI / 72.0);
//					    paper.setImageableArea(0, 0, paper.getWidth(), paper.getHeight());
//					    PageFormat pf = new PageFormat();
//					    pf.setOrientation(orientation);
//					    pf.setPaper(paper);
//
//						BufferedImage img = new BufferedImage((int)pf.getWidth(), (int)pf.getHeight(), BufferedImage.TYPE_INT_RGB);
//						final Graphics2D g2d = img.createGraphics();
//
////						ProxyFactory factory = new ProxyFactory();
////					    factory.setSuperclass(Graphics2D.class);
////					    factory.setFilter(
////					            new MethodFilter() {
////					                @Override
////					                public boolean isHandled(Method method) {
////					                    return Modifier.isAbstract(method.getModifiers());
////					                }
////					            }
////					        );
////
////					    MethodHandler handler = new MethodHandler() {
////					        @Override
////					        public Object invoke(Object self, Method thisMethod, Method proceed, Object[] args) throws Throwable {
//////					            System.out.println("Handling " + thisMethod + " via the method handler");
////								if (thisMethod.getName() == "setFont") {
////									customSetFont((Font)args[0]);
////									return null;
////								}
////								return thisMethod.invoke(g2d, args);
////					        }
////							public void customSetFont(Font font) {
////								float fontSize = font.getSize2D();
////								Font newFont = font.deriveFont(fontSize * DPI / 72.0f);
////								g2d.setFont(newFont);
////							}
////					    };
////
////					    Graphics2D customG2d;
////						try {
////							customG2d = (Graphics2D)factory.create(new Class<?>[0], new Object[0], handler);
////			                chartPanels.get(0).print(customG2d, pf, 0);
////							g2d.dispose();
////						} catch (NoSuchMethodException e2) {
////							// TODO Auto-generated catch block
////							e2.printStackTrace();
////						} catch (IllegalArgumentException e2) {
////							// TODO Auto-generated catch block
////							e2.printStackTrace();
////						} catch (InstantiationException e2) {
////							// TODO Auto-generated catch block
////							e2.printStackTrace();
////						} catch (IllegalAccessException e2) {
////							// TODO Auto-generated catch block
////							e2.printStackTrace();
////						} catch (InvocationTargetException e2) {
////							// TODO Auto-generated catch block
////							e2.printStackTrace();
////						}
//
////						InvocationHandler handler = new InvocationHandler() {
////							@Override
////							public Object invoke(Object proxy, Method method, Object[] args)
////									throws Throwable {
////								if (method.getName() == "setFont") {
////									setFont((Font)args[0]);
////									return null;
////								}
////								return method.invoke(g2d, args);
////							}
////							public void setFont(Font font) {
////								float fontSize = font.getSize2D();
////								Font newFont = font.deriveFont(fontSize * DPI / 72.0f);
////								g2d.setFont(newFont);
////							}
////						};
////						Graphics2D customG2d = (Graphics2D)Proxy.newProxyInstance(
////								Graphics2D.class.getClassLoader(),
////                                new Class[] { Graphics2D.class },
////                                handler);
////						Graphics2D customG2d = new  {
////							@Override
////							public void setFont(Font font) {
////								float fontSize = font.getSize2D();
////								Font newFont = font.deriveFont(fontSize * DPI / 72.0f);
////								super.setFont(newFont);
////							}
////						};
////		                chartPanels.get(0).print(customG2d, pf, 0);
////						g2d.dispose();
//
////						File outputFile = new File("/home/bhepp/abc.png");
////						try {
////							outputFile.createNewFile();
////							ImageIO.write(img, "png", outputFile);
////						} catch (IOException e1) {
////							String errorMsg = "Failed to save plots as PNG:\n" + e1.toString();
////							JOptionPane.showMessageDialog(GridWindow.this, errorMsg, "Error", ERROR);
////						}
//
////					    int     res = pj.getPageResolution();
////					      Dimension d = pj.getPageDimension();
////					      System.out.println( "Resolution : " + res + "\n" +
////					                          "Width : " + d.width + "\n" +
////					                          "Height : " + d.height + "\n" +
////					                          "Pixel on page : " + (res * d.width * d.height) );
//
////					    pj
////					    GridWindow.this.print(g, pj., page)
//
////						DocFlavor flavor = DocFlavor.INPUT_STREAM.PDF;
////						aset.add(MediaSizeName.ISO_A4);
////						PrintService[] service = PrintServiceLookup.lookupPrintServices(flavor, aset);
////						PrinterJob job = PrinterJob.getPrinterJob();
////						job.setPrintService(service[0]);
////						PrinterJob job = PrinterJob.getPrinterJob();
////						job.setPrintable(GridWindow.this);
////						boolean doPrint = job.printDialog(aset);
////						PageFormat pf = job.getPageFormat(aset);
////						pf = pf;
////						if (doPrint) {
////							try {
////								job.print(aset);
////							} catch (PrinterException e1) {
////								String errorMsg = "Failed to print plots:\n" + e1.toString();
////								JOptionPane.showMessageDialog(GridWindow.this, errorMsg, "Error", ERROR);
////							}
////						}
//					}
//				}
//			);
//		// Print as PDF entry
//		JMenuItem saveAsPDFEntry = new JMenuItem("Print as PDF");
//		saveAsPDFEntry.setAccelerator(KeyStroke.getKeyStroke('P', InputEvent.CTRL_MASK));
//		saveAsPDFEntry.addActionListener(
//				new ActionListener() {
//					@Override
//					public void actionPerformed(ActionEvent e) {
////						JFileChooser fc = new CustomJFileChooser();
////						int returnValue = fc.showSaveDialog(GridWindow.this);
////						if (returnValue == JFileChooser.APPROVE_OPTION) {
////							File outputFile = fc.getSelectedFile();
////							try {
////								outputFile.createNewFile();
////								PDFJob job = new PDFJob(new FileOutputStream(outputFile));
////								Graphics2D g2d = (Graphics2D)job.getGraphics(PageFormat.PORTRAIT);
////							    Dimension d = job.getPageDimension();
////							    Rectangle area = job.getCurrentPage().getImageableArea();
////								drawPlots(g2d, area.getX(), area.getY(), area.getWidth(), area.getHeight());
////								g2d.dispose();
////								job.end();
////							} catch (IOException e1) {
////								String errorMsg = "Failed to save plots as PDF:\n" + e1.toString();
////								JOptionPane.showMessageDialog(GridWindow.this, errorMsg, "Error", ERROR);
////							}
////						}
//						JFileChooser fc = new CustomJFileChooser();
//						fc.setFileFilter(new FileNameExtensionFilter("PDF documents", "pdf"));
//						int returnValue = fc.showSaveDialog(GridWindow.this);
//						if (returnValue == JFileChooser.APPROVE_OPTION) {
//							File outputFile = fc.getSelectedFile();
//							Document document = null;
//							PdfWriter writer = null;
//							try {
//								document = new Document(PageSize.A4.rotate());
//								outputFile.createNewFile();
//								writer = PdfWriter.getInstance(document, new FileOutputStream(outputFile));
//								document.open();
//								PdfContentByte canvas = writer.getDirectContent();
//
////								// These lines do not influence the pdf document,
////								// but are there to tell the Printable how to print
//				                final int DPI = 300;
//								final double a4WidthInch = 11.6929134; // Equals 210mm
//				                final double a4HeightInch = 8.26771654; // Equals 297mm
//				                int orientation = PageFormat.PORTRAIT;
//				                Paper paper = new Paper();
//				                paper.setSize(a4WidthInch * DPI, a4HeightInch * DPI);
//				                paper.setImageableArea(
//				                		DPI, DPI,
//				                		a4WidthInch * DPI - 2 * DPI,
//				                		a4HeightInch * DPI - 2 * DPI); // 1 inch margins
//				                PageFormat pageFormat = new PageFormat();
//				                pageFormat.setPaper(paper);
//				                pageFormat.setOrientation(orientation);
////
////				                double width = (pageFormat.getWidth());
////				                double height = (pageFormat.getHeight());
////				                canvas.saveState();
////				                canvas.setColorStroke(GrayColor.BLACK);
////				                canvas.setColorFill(GrayColor.BLACK);
////				                canvas.circle(200, 200, 50);
////								int x0 = (int)paper.getImageableX();
////								int y0 = (int)paper.getImageableY();
////								int w = (int)paper.getImageableWidth();
////								int h = (int)paper.getImageableHeight();
//								com.itextpdf.text.Rectangle pageSize = document.getPageSize();
//								double scaleX = 72.0 / DPI;
//								double scaleY = 72.0 / DPI;
//								double graphicsWidth = pageSize.getWidth() / scaleX;
//								double graphicsHeight = pageSize.getHeight() / scaleY;
//
//								BufferedImage img = new BufferedImage((int)graphicsWidth, (int)graphicsHeight, BufferedImage.TYPE_INT_RGB);
//								Graphics2D g2d = img.createGraphics();
//								drawPlots(g2d, 0, 0, img.getWidth(), img.getHeight());
//
//				                g2d = new PdfGraphics2D(
//				                		canvas, pageSize.getWidth(), pageSize.getHeight(), new PdfFontMapper());
////				                g2d = new PdfGraphics2D(
////				                		canvas, (float)graphicsWidth, (float)graphicsHeight, new PdfFontMapper());
////				                g2d.drawImage(img, new AffineTransform(scaleX, 0, 0, scaleY, 0, 0), null);
//				                g2d.scale(scaleX, scaleY);
//				                chartPanels.get(0).print(g2d, pageFormat, 0);
////								GridWindow.this.drawPlots(g2d, 0, 0, (int)graphicsWidth, (int)graphicsHeight);
//								g2d.dispose();
////								canvas.fillStroke();
////				                canvas.restoreState();
//
//							} catch (DocumentException e1) {
//								String errorMsg = "Failed to save plots as PDF:\n" + e1.toString();
//								JOptionPane.showMessageDialog(GridWindow.this, errorMsg, "Error", ERROR);
//							} catch (IOException e1) {
//								String errorMsg = "Failed to save plots as PDF:\n" + e1.toString();
//								JOptionPane.showMessageDialog(GridWindow.this, errorMsg, "Error", ERROR);
//							} finally {
//								if (document != null)
//									document.close();
//								if (writer != null)
//									writer.close();
//							}
//						}
//					}
//				}
//			);
//		// Save as entry
//		JMenuItem saveAsPNGEntry = new JMenuItem("Save as PNG...");
//		saveAsPNGEntry.setAccelerator(KeyStroke.getKeyStroke('S', InputEvent.CTRL_MASK));
//		saveAsPNGEntry.addActionListener(
//				new ActionListener() {
//					@Override
//					public void actionPerformed(ActionEvent e) {
//						JFileChooser fc = new CustomJFileChooser();
//						int returnValue = fc.showSaveDialog(GridWindow.this);
//						if (returnValue == JFileChooser.APPROVE_OPTION) {
//							final int width = GridWindow.this.cols * GridWindow.this.getSavePlotWidth()
//									+ (GridWindow.this.cols - 1) * GridWindow.this.getSaveHorizontalSpacing();
//							final int height = GridWindow.this.rows * GridWindow.this.getSavePlotHeight()
//									+ (GridWindow.this.rows - 1) * GridWindow.this.getSaveVerticalSpacing();
//							BufferedImage outputImg = new BufferedImage(width, height, BufferedImage.TYPE_INT_BGR);
//							Graphics2D g2d = outputImg.createGraphics();
//							drawPlots(g2d, 0, 0, outputImg.getWidth(), outputImg.getHeight());
////							g2d.setBackground(Color.lightGray);
////							g2d.clearRect(0, 0, outputImg.getWidth(), outputImg.getHeight());
////							int row = 0;
////							int col = 0;
////							for (ChartPanel panel : GridWindow.this.chartPanels) {
////								JFreeChart chart = panel.getChart();
////								int x = col * (GridWindow.this.getSavePlotWidth() + GridWindow.this.getSaveHorizontalSpacing());
////								int y = row * (GridWindow.this.getSavePlotHeight() + GridWindow.this.getSaveVerticalSpacing());
////								Rectangle area = new Rectangle(x, y, GridWindow.this.getSavePlotWidth(), GridWindow.this.getSavePlotHeight());
////								chart.draw(g2d, area);
////								col++;
////								if (col >= cols) {
////									col = 0;
////									row++;
////								}
////							}
//							File outputFile = fc.getSelectedFile();
//							try {
//								ImageIO.write(outputImg, "png", outputFile);
//							} catch (IOException e1) {
//								String errorMsg = "Failed to save plots as PNG:\n" + e1.toString();
//								JOptionPane.showMessageDialog(GridWindow.this, errorMsg, "Error", ERROR);
//							}
//						}
//					}
//				}
//			);
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
//		fileMenu.add(saveAsPdfWithREntry);
		fileMenu.add(exportToMatlabFileEntry);
		fileMenu.add(plotInMatlabEntry);
//		fileMenu.add(testEntry);
//		fileMenu.add(printEntry);
//		fileMenu.add(saveAsPDFEntry);
//		fileMenu.add(saveAsPNGEntry);
		fileMenu.addSeparator();
		fileMenu.add(quitEntry);
		menubar.add(fileMenu);
		return menubar;
	}

//	public void drawPlots(Graphics2D g2d, int x0, int y0, int width, int height) {
//		g2d.setBackground(Color.white);
//		g2d.clearRect(x0, y0, width, height);
//		int plotWidth = (int)Math.round(width / (double)cols);
//		int plotHeight = (int)Math.round(height / (double)rows);
//		int row = 0;
//		int col = 0;
//		for (ChartPanel panel : GridWindow.this.chartPanels) {
//			JFreeChart chart = panel.getChart();
//			int x = x0 + col * plotWidth;
//			int y = y0 + row * plotHeight;
//			java.awt.Rectangle area = new java.awt.Rectangle(x, y, plotWidth, plotHeight);
//			chart.draw(g2d, area);
//			col++;
//			if (col >= cols) {
//				col = 0;
//				row++;
//			}
//		}
//	}

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
