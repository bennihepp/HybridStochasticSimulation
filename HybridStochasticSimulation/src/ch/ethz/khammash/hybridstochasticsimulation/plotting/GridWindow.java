package ch.ethz.khammash.hybridstochasticsimulation.plotting;

import static com.google.common.base.Preconditions.*;

import java.awt.Component;
import java.awt.GridLayout;

import javax.swing.JFrame;

import org.jfree.ui.ApplicationFrame;

public class GridWindow extends ApplicationFrame {
	private static final long serialVersionUID = 6566839328303930162L;

	public GridWindow(String title, int rows, int cols) {
		super(title);
		checkArgument(rows > 0);
		checkArgument(cols > 0);
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
