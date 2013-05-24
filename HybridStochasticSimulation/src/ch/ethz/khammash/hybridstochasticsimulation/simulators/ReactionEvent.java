package ch.ethz.khammash.hybridstochasticsimulation.simulators;

import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.RealVector;

public class ReactionEvent implements Comparable<Double> {

	private int reaction;
	private double time;
//	private double[] newX;
	private RealVector newX;

	public ReactionEvent(double t) {
//		this(-1, t, new double[0]);
		this(-1, t, new ArrayRealVector());
	}

	public ReactionEvent(int reaction, double t, double[] newX) {
		this.reaction = reaction;
		this.time = t;
		this.newX = new ArrayRealVector(newX);
	}

	public ReactionEvent(int reaction, double t, RealVector newX) {
		this.reaction = reaction;
		this.time = t;
		this.newX = newX.copy();
	}

	public int getReaction() {
		return reaction;
	}

	public double getTime() {
		return time;
	}

	public double[] getNewX() {
		return newX.toArray();
	}

	public RealVector getNewXVector() {
		return newX;
	}

	@Override
	public int compareTo(Double t) {
		return Double.compare(time, t);
	}

}