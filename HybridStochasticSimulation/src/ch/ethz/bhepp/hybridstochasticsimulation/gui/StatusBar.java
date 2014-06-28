package ch.ethz.bhepp.hybridstochasticsimulation.gui;

import java.awt.Dimension;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.util.Deque;
import java.util.LinkedList;

import javax.swing.BoxLayout;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.SwingConstants;
import javax.swing.Timer;
import javax.swing.border.BevelBorder;

public class StatusBar extends JPanel {

	private static final long serialVersionUID = -8534232086639359456L;

	private JLabel statusLabel;
	private Deque<TemporaryStatusText> statusTextDeque;
	private String permanentText;

	public StatusBar() {
		setBorder(new BevelBorder(BevelBorder.LOWERED));
		Dimension pf = getPreferredSize();
		pf.setSize(pf.getWidth(), 16.0);
		setPreferredSize(pf);
		setLayout(new BoxLayout(this, BoxLayout.X_AXIS));
		statusLabel = new JLabel();
		statusLabel.setHorizontalAlignment(SwingConstants.LEFT);
		add(statusLabel);
		statusTextDeque = new LinkedList<TemporaryStatusText>();
	}

	public void setText(String text) {
		synchronized (statusTextDeque) {
			statusLabel.setText(text);
			permanentText = text;
		}
	}

	public void clearText() {
		setText("");
	}

	public void setText(String text, int duration) {
		setText(text, duration, false);
	}

	public void setText(String text, int duration, boolean clearAfterwards) {
		synchronized (statusTextDeque) {
			if (statusTextDeque.isEmpty()) {
				permanentText = statusLabel.getText();
				statusLabel.setText(text);
				ActionListener listener = new ActionListener() {
					@Override
					public void actionPerformed(ActionEvent e) {
						synchronized (statusTextDeque) {
							statusTextDeque.removeLast();
							if (statusTextDeque.size() > 0) {
								TemporaryStatusText tst = statusTextDeque.getLast();
								statusLabel.setText(tst.getText());
								startTimerWithoutRepeats(tst.getDuration(), this);
							} else
								statusLabel.setText(permanentText);
						}
					}
				};
				startTimerWithoutRepeats(duration, listener);
			}
			if (clearAfterwards)
				permanentText = "";
			statusTextDeque.addFirst(new TemporaryStatusText(text, duration));
		}
	}

	private void startTimerWithoutRepeats(int duration, ActionListener listener) {
		Timer t = new Timer(duration, listener);
		t.setRepeats(false);
		t.start();
	}

	public class TemporaryStatusText {

		private String text;
		private int duration;

		public TemporaryStatusText(String text, int duration) {
			this.text = text;
			this.duration = duration;
		}

		public String getText() {
			return text;
		}

		public int getDuration() {
			return duration;
		}

	}

}
