package ch.ethz.bhepp.hybridstochasticsimulation.outputs;

import static com.google.common.base.Preconditions.checkArgument;

import java.util.ArrayList;
import java.util.List;

import ch.ethz.bhepp.hybridstochasticsimulation.trajectories.FiniteTrajectory;

public class FixedTimepointOutput implements OutputRecord {

	private static final String T_RECORD_LABEL = "t";
	private static final String X_RECORD_LABEL = "x";

	private static class FixedTimepointState {
		private String label;
		private double[] record;
	}

	private static class FixedTimepointOutputChild implements OutputRecord {

		private List<FixedTimepointState> recordList;

		public FixedTimepointOutputChild(double time, FiniteTrajectory tr, List<String> labels) {
			checkArgument(tr.getNumberOfStates() == labels.size());
			double[] state = tr.getInterpolatedState(time);
			recordList = new ArrayList<>(tr.getNumberOfStates());
			for (int i=0; i < tr.getNumberOfStates(); i++) {
				FixedTimepointState rec = new FixedTimepointState();
				rec.label = labels.get(i);
				rec.record = new double[1];
				rec.record[0] = state[i];
				recordList.add(rec);
			}
		}

		@Override
		public int getNumberOfRecords() {
			return recordList.size();
		}

		@Override
		public String getRecordLabel(int recordIndex) {
			return recordList.get(recordIndex).label;
		}

		@Override
		public double[] getRecord(int recordIndex) {
			return recordList.get(recordIndex).record;
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

	private FixedTimepointOutputChild child;
	private double[] timeArray;

	public FixedTimepointOutput(double time, FiniteTrajectory tr, List<String> labels) {
		checkArgument(tr.getNumberOfStates() == labels.size());
		this.child = new FixedTimepointOutputChild(time, tr, labels);
		this.timeArray = new double[1];
		this.timeArray[0] = time;
	}

	@Override
	public int getNumberOfRecords() {
		return 1;
	}

	@Override
	public String getRecordLabel(int recordIndex) {
		if (recordIndex == 0)
			return T_RECORD_LABEL;
		else
			throw new IndexOutOfBoundsException();
	}

	@Override
	public double[] getRecord(int recordIndex) {
		if (recordIndex == 0)
			return timeArray;
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
			return X_RECORD_LABEL;
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
