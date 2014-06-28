package ch.ethz.bhepp.hybridstochasticsimulation.simulators.ode;

import java.util.List;

import ch.ethz.bhepp.hybridstochasticsimulation.simulators.PDMPEventObserver;
import ch.ethz.bhepp.ode.EventFunction;
import ch.ethz.bhepp.ode.EventObserver;

public class EventObserverAdapter implements EventFunction, EventObserver {

	private PDMPEventObserver[] eventObservers;
	private long eventCount;

	public EventObserverAdapter(List<PDMPEventObserver> eventObservers) {
		this.eventObservers = eventObservers.toArray(new PDMPEventObserver[0]);
	}

	@Override
	public int getNumberOfEventValues() {
		return eventObservers.length;
	}

	@Override
	public void computeEventValues(double t, double[] x, double[] values) {
		for (int i=0; i < eventObservers.length; i++)
			values[i] = eventObservers[i].g(t, x);
	}

	@Override
	public void initialize(double t0, double[] x0, double t1) {
		eventCount = 0;
	}

	@Override
	public EventAction report(int eventIndex, double t, double[] x) {
		eventCount++;
		eventObservers[eventIndex].eventOccurred(t, x, false);
		return EventAction.STOP;
	}

	public long getEventCount() {
		return eventCount;
	}
}
