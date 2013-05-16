package ch.ethz.khammash.hybridstochasticsimulation;

import java.awt.Component;
import java.awt.GridLayout;

import javax.swing.JFrame;

import org.jfree.ui.ApplicationFrame;

public class GridWindow extends ApplicationFrame {
	private static final long serialVersionUID = 6566839328303930162L;

	public GridWindow(String title, int rows, int cols) {
		super(title);
		setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		setLayout(new GridLayout(rows, cols));
	}

	@Override
	public Component add(Component comp) {
		Component q = super.add(comp);
		pack();
		return q;
	}
}
