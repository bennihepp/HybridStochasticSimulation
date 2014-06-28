package ch.ethz.bhepp.hybridstochasticsimulation.simulators;

import ch.ethz.bhepp.hybridstochasticsimulation.models.ReactionNetworkModel;

public abstract class AbstractSimulator<T extends ReactionNetworkModel> implements Simulator<T> {

	boolean printMessages;
	boolean showProgress;

	public void setPrintMessages(boolean printMessages) {
		this.printMessages = printMessages;
	}

	public void setShowProgress(boolean showProgress) {
		this.showProgress = showProgress;
	}

}
