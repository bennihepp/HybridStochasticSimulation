package ch.ethz.bhepp.hybridstochasticsimulation.simulators;

import java.util.LinkedList;
import java.util.List;

import ch.ethz.bhepp.hybridstochasticsimulation.models.StateBoundEventListener;

public class PDMPEventObserverCollector implements PDMPEventObserver {

	private StateBoundEventListener listener;
	private List<PDMPEventObserver> eventObservers;

	public PDMPEventObserverCollector(StateBoundEventListener listener) {
		this.listener = listener;
		this.eventObservers = new LinkedList<PDMPEventObserver>();
	}

	public PDMPEventObserverCollector(StateBoundEventListener listener, List<PDMPEventObserver> eventObservers) {
		this.listener = listener;
		this.eventObservers = eventObservers;
	}

	public void add(PDMPEventObserver observer) {
		eventObservers.add(observer);
	}

	public void remove(PDMPEventObserver observer) {
		eventObservers.remove(observer);
	}

	@Override
	public void init(double t0, double[] x0, double t) {
	}

	@Override
	public double g(double t, double[] x) {
		double value = 1.0;
		for (PDMPEventObserver observer : eventObservers)
			value *= observer.g(t, x);
		return value;
	}

	@Override
	public Action eventOccurred(double t, double[] x, boolean increasing) {
		listener.stateBoundEventOccured(t, x);
		return Action.STOP;
	}

	@Override
	public void resetState(double t, double[] x) {
	}

}
