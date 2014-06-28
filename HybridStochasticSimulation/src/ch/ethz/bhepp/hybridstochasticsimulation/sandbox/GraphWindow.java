package ch.ethz.khammash.hybridstochasticsimulation.sandbox;


import java.awt.Color;
import java.awt.Dimension;
import java.awt.Rectangle;
import java.awt.geom.Rectangle2D;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Set;

import javax.swing.JFrame;
import javax.swing.JScrollPane;

import org.apache.commons.math3.util.FastMath;
import org.jgraph.JGraph;
import org.jgraph.graph.AttributeMap;
import org.jgraph.graph.DefaultGraphCell;
import org.jgraph.graph.GraphConstants;
import org.jgrapht.ext.JGraphModelAdapter;

import ch.ethz.khammash.hybridstochasticsimulation.graphs.ReactionEdge;
import ch.ethz.khammash.hybridstochasticsimulation.graphs.DefaultReactionNetworkGraph;
import ch.ethz.khammash.hybridstochasticsimulation.graphs.SpeciesVertex;


/**
 * Adapted from jgrapht example from Barak Naveh, Aug 3, 2003.
 */

public class GraphWindow extends JFrame {

	private static final long serialVersionUID = -8373894585699085082L;

	private static final Color     DEFAULT_BG_COLOR = Color.decode( "#FAFBFF" );
    private static final Dimension DEFAULT_SIZE = new Dimension( 530, 320 );

    private JScrollPane scrollPane;
    private DefaultReactionNetworkGraph graph;
    private JGraphModelAdapter<SpeciesVertex,ReactionEdge> jgAdapter;
    private JGraph jgraph;

    /**
     * @see java.applet.Applet#init().
     */
    public GraphWindow(DefaultReactionNetworkGraph graph) {
    	setPreferredSize(new Dimension(1024, 786));

    	this.graph = graph;

        // Create a visualization using JGraph, via an adapter
        jgAdapter = new JGraphModelAdapter<SpeciesVertex,ReactionEdge>(graph);
        jgraph = new JGraph(jgAdapter);
        adjustDisplaySettings(jgraph);

        scrollPane = new JScrollPane(jgraph);
        getContentPane().add(scrollPane);
    }

//    public GraphWindow(DirectedGraph<Integer,DefaultEdge> graph) {
//    	this.graph = graph;
//
//        // Create a visualization using JGraph, via an adapter
//        jgAdapter = new JGraphModelAdapter<Integer,DefaultEdge>(graph);
//        jgraph = new JGraph(jgAdapter);
//        adjustDisplaySettings(jgraph);
//
//        scrollPane = new JScrollPane(jgraph);
//        getContentPane().add(scrollPane);
//    }

    public void autoPosition() {
//    	Object[] roots = new Object[0];
//    	JGraphFacade facade = new JGraphFacade(jgAdapter, roots, false, false, false, true, new JGraphConstantCostFunction(1.0), JGraphAlgebra.getSharedInstance());
//    	facade.circle(facade.getVertices());
//    	JGraphOrganicLayout graphLayout = new JGraphOrganicLayout();
//    	graphLayout.run(new JGraphFacade(jgraph));
//        jgraph = new JGraph(jgAdapter);
//        scrollPane.setViewportView(jgraph);
        // position vertices nicely within JGraph component
        Dimension size = getSize();
        Set<SpeciesVertex> vertexSet = graph.vertexSet();
        Iterator<SpeciesVertex> it = vertexSet.iterator();
        for (int i=0; i < vertexSet.size(); i++) {
        	double phi = 2 * FastMath.PI * i / (double)vertexSet.size();
        	double x = size.getWidth() / 3.0 * FastMath.cos(phi) + size.getWidth() / 2.0;
        	double y = size.getHeight() / 3.0 * FastMath.sin(phi) + size.getHeight() / 2.0;
        	positionVertexAt(it.next(), x, y);
        }
        revalidate();
    }


    private void adjustDisplaySettings(JGraph jg) {
        jg.setPreferredSize(DEFAULT_SIZE);
        Color c = DEFAULT_BG_COLOR;
        jg.setBackground(c);
    }

    public void positionVertexAt(SpeciesVertex vertex, double x, double y) {
        DefaultGraphCell cell = jgAdapter.getVertexCell(vertex);
        AttributeMap attr = cell.getAttributes();
        Rectangle2D b = GraphConstants.getBounds(attr);

        GraphConstants.setBounds(attr, new Rectangle((int)x, (int)y, (int)b.getWidth(), (int)b.getHeight()));

        HashMap<DefaultGraphCell,AttributeMap> cellAttr = new HashMap<DefaultGraphCell,AttributeMap>();
        cellAttr.put(cell, attr);
        jgAdapter.edit(cellAttr, null, null, null);
    }

//    public static void main(String[] args) {
//        // create a JGraphT graph
//    	DirectedGraph<Integer,DefaultEdge> graph = new DefaultDirectedGraph<Integer,DefaultEdge>(DefaultEdge.class);
//
//        // add some sample data (graph manipulated via JGraphT)
//        graph.addVertex(0);
//        graph.addVertex(1);
//        graph.addVertex(2);
//        graph.addVertex(3);
//
//        graph.addEdge(0, 1);
//        graph.addEdge(1, 2);
//        graph.addEdge(2, 3);
//        graph.addEdge(3, 0);
//
//		GraphWindow window = new GraphWindow(graph);
//		window.pack();
//		window.setVisible(true);
//		window.autoPosition();
//    }

}
