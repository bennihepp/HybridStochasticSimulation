package ch.ethz.khammash.hybridstochasticsimulation.outputs;


public interface OutputRecord {

	int getNumberOfRecords();

	String getRecordLabel(int record);

	double[] getRecord(int record);

	int getNumberOfChildren();

	String getChildLabel(int child);

	OutputRecord getChild(int child);

}
