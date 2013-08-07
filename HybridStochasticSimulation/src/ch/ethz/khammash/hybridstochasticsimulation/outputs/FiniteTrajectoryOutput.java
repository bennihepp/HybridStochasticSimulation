package ch.ethz.khammash.hybridstochasticsimulation.outputs;

import static com.google.common.base.Preconditions.checkArgument;

import java.util.List;

import ch.ethz.khammash.hybridstochasticsimulation.trajectories.FiniteTrajectory;

public class FiniteTrajectoryOutput implements OutputRecord {

	private static final String TSERIES_RECORD_LABEL = "tSeries";
	private static final String XSERIES_RECORD_LABEL = "xSeries";

	private static class FiniteTrajectoryOutputChild implements OutputRecord {

		private FiniteTrajectory tr;
		private List<String> labels;

		public FiniteTrajectoryOutputChild(FiniteTrajectory tr, List<String> labels) {
			checkArgument(tr.getNumberOfStates() == labels.size());
			this.tr = tr;
			this.labels = labels;
		}

		@Override
		public int getNumberOfRecords() {
			return tr.getNumberOfStates();
		}

		@Override
		public String getRecordLabel(int recordIndex) {
			return labels.get(recordIndex);
		}

		@Override
		public double[] getRecord(int recordIndex) {
			return tr.getxSeries(recordIndex);
		}

		@Override
		public int getNumberOfChildren() {
			return 0;
		}

		@Override
		public String getChildLabel(int childIndex) {
			throw new IndexOutOfBoundsException();
		}

		@Override
		public OutputRecord getChild(int child) {
			throw new IndexOutOfBoundsException();
		}

	}

	private FiniteTrajectoryOutputChild child;
	private double[] tSeries;

	public FiniteTrajectoryOutput(FiniteTrajectory tr, List<String> labels) {
		checkArgument(tr.getNumberOfStates() == labels.size());
		this.child = new FiniteTrajectoryOutputChild(tr, labels);
		tSeries = tr.gettSeries();
	}

	@Override
	public int getNumberOfRecords() {
		return 1;
	}

	@Override
	public String getRecordLabel(int recordIndex) {
		if (recordIndex == 0)
			return TSERIES_RECORD_LABEL;
		else
			throw new IndexOutOfBoundsException();
	}

	@Override
	public double[] getRecord(int recordIndex) {
		if (recordIndex == 0)
			return tSeries;
		else
			throw new IndexOutOfBoundsException();
	}

	@Override
	public int getNumberOfChildren() {
		return 1;
	}

	@Override
	public String getChildLabel(int childIndex) {
		if (childIndex == 0)
			return XSERIES_RECORD_LABEL;
		else
			throw new IndexOutOfBoundsException();
	}

	@Override
	public OutputRecord getChild(int childIndex) {
		if (childIndex == 0)
			return child;
		else
			throw new IndexOutOfBoundsException();
	}

}
