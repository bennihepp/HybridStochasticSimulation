/**
 * @author      Benjamin Hepp <benjamin.hepp@bsse.ethz.ch>
 */
package ch.ethz.bhepp.ode;


public class FiniteTimepointProvider implements TimepointProvider {

	private double[] tSeries;
	private int nextIndex;
	private double tCurrent;

	public FiniteTimepointProvider(double t0, double t1) {
		tSeries = new double[2];
		tSeries[0] = t0;
		tSeries[1] = t1;
		tCurrent = t0;
		nextIndex = 0;
	}

	public FiniteTimepointProvider(double[] tSeries) {
		this.tSeries = tSeries;
		tCurrent = tSeries[0];
		nextIndex = 0;
	}

	@Override
	public double getInitialTimepoint() {
		return tSeries[0];
	}

	@Override
	public double getLastTimepoint() {
		return tSeries[tSeries.length - 1];
	}

	@Override
	public void reset() {
		nextIndex = 0;
	}

	public void setCurrentTimepoint(double t) {
		this.tCurrent = t;
	}

	@Override
	public double getCurrentTimepoint() {
		return tCurrent;
	}

	@Override
	public boolean hasNextTimepoint(double tCurrent) {
		while (nextIndex < tSeries.length)
			if (tCurrent >= tSeries[nextIndex])
				nextIndex++;
			else
				return true;
		return false;
//		if (nextIndex >= tSeries.length)
//			return false;
//		double tNext = peekNextTimepoint();
//		return tCurrent < tNext;
	}

	@Override
	public double getNextTimepoint(double tCurrent) throws NoMoreTimepointsException {
		if (hasNextTimepoint(tCurrent)) {
			this.tCurrent = tCurrent;
			return tSeries[nextIndex];
		} else
			throw new NoMoreTimepointsException("No more timepoints available");
//		if (tCurrent >= tSeries[nextIndex])
//			if (nextIndex + 1 < tSeries.length)
//		if (nextIndex < tSeries.length) {
//			double t = tSeries[nextIndex];
//			nextIndex++;
//			return t;
//		}
//		throw new NoMoreTimepointsException("No more timepoints available");
	}

//	private double peekNextTimepoint() {
//		return tSeries[nextIndex];
//	}

}
