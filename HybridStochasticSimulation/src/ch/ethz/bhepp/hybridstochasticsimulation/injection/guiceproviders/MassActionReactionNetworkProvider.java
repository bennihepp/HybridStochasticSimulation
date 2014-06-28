package ch.ethz.bhepp.hybridstochasticsimulation.injection.guiceproviders;

import java.io.File;
import java.io.IOException;

import javax.xml.parsers.ParserConfigurationException;

import org.apache.commons.configuration.HierarchicalConfiguration;
import org.xml.sax.SAXException;

import ch.ethz.bhepp.hybridstochasticsimulation.examples.ExampleConfigurationFactory;
import ch.ethz.bhepp.hybridstochasticsimulation.io.StochKitNetworkReader;
import ch.ethz.bhepp.hybridstochasticsimulation.io.StochKitNetworkReader.FileFormatException;
import ch.ethz.bhepp.hybridstochasticsimulation.networks.MassActionReactionNetwork;

import com.google.inject.Inject;

public class MassActionReactionNetworkProvider extends AbstractProvider<MassActionReactionNetwork> {

	@Inject
	public MassActionReactionNetworkProvider(HierarchicalConfiguration config) {
		super(config, "ModelParameters");
	}

	@Override
	public MassActionReactionNetwork get() {
		String inputFileName = config().getString("inputFile", null);
		if (inputFileName != null) {
			File inputFile = new File(inputFileName);
			try {
				MassActionReactionNetwork network = StochKitNetworkReader.readUnaryBinaryNetworkFromFile(inputFile);
				return network;
			} catch (ParserConfigurationException | SAXException | IOException | FileFormatException e) {
				throw new RuntimeException(e);
			}
		}
		String exampleName = config().getString("example");
		return ExampleConfigurationFactory.getInstance().createExampleConfiguration(exampleName).net;
	}

}
