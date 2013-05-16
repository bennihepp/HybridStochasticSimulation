public class ReactionRecord implements Comparable<Double> {

	protected int reaction;
	protected double time;
	protected double[] newX;

	public ReactionRecord(double t) {
		this(-1, t, new double[0]);
	}

	public ReactionRecord(int reaction, double t, double[] newX) {
		this.reaction = reaction;
		this.time = t;
		this.newX = newX;
	}

	public int getReaction() {
		return reaction;
	}

	public double getTime() {
		return time;
	}

	public double[] getNewX() {
		return newX;
	}

	@Override
	public int compareTo(Double t) {
		return Double.compare(time, t);
	}

}