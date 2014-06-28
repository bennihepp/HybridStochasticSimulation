package ch.ethz.bhepp.hybridstochasticsimulation.matlab;

import matlabcontrol.MatlabConnectionException;
import matlabcontrol.MatlabInvocationException;
import matlabcontrol.MatlabProxy;
import matlabcontrol.MatlabProxyFactory;
import matlabcontrol.MatlabProxyFactoryOptions;
import matlabcontrol.MatlabProxyFactoryOptions.Builder;

public class MatlabSession {

	MatlabProxy proxy;
	boolean hidden;

	public MatlabSession() {
		this(false);
	}

	public MatlabSession(boolean hidden) {
		this.hidden = hidden;
	}

	public void start() throws MatlabConnectionException {
		// Create a proxy, which we will use to control MATLAB
		Builder builder = new MatlabProxyFactoryOptions.Builder();
		builder.setHidden(hidden);
		builder.setUsePreviouslyControlledSession(true);
		MatlabProxyFactoryOptions options = builder.build();
		MatlabProxyFactory factory = new MatlabProxyFactory(options);
		proxy = factory.getProxy();
	}

	public void stop() throws MatlabInvocationException {
		if (isExistingSession())
		    // Disconnect the proxy from MATLAB
		    proxy.disconnect();
		else
		    // Exit MATLAB session
		    proxy.exit();
	}

	public MatlabProxy getProxy() {
		return proxy;
	}

	public boolean isExistingSession() {
		return proxy.isExistingSession();
	}

}
