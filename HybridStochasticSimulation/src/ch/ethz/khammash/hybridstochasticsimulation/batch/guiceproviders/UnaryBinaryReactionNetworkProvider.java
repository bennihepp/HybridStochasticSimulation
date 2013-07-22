package ch.ethz.khammash.hybridstochasticsimulation.batch.guiceproviders;

import java.io.File;
import java.io.IOException;

import javax.xml.parsers.ParserConfigurationException;

import org.apache.commons.configuration.HierarchicalConfiguration;
import org.xml.sax.SAXException;

import ch.ethz.khammash.hybridstochasticsimulation.examples.ExampleConfigurationFactory;
import ch.ethz.khammash.hybridstochasticsimulation.io.StochKitNetworkReader;
import ch.ethz.khammash.hybridstochasticsimulation.io.StochKitNetworkReader.FileFormatException;
import ch.ethz.khammash.hybridstochasticsimulation.networks.UnaryBinaryReactionNetwork;

import com.google.inject.Inject;

public class UnaryBinaryReactionNetworkProvider extends AbstractProvider<UnaryBinaryReactionNetwork> {

	@Inject
	public UnaryBinaryReactionNetworkProvider(HierarchicalConfiguration config) {
		super(config, "ModelParameters");
	}

	@Override
	public UnaryBinaryReactionNetwork get() {
		String inputFileName = config().getString("inputFile", null);
		if (inputFileName != null) {
			File inputFile = new File(inputFileName);
			try {
				UnaryBinaryReactionNetwork network = StochKitNetworkReader.readUnaryBinaryNetworkFromFile(inputFile);
				return network;
			} catch (ParserConfigurationException | SAXException | IOException | FileFormatException e) {
				throw new RuntimeException(e);
			}
		}
		String exampleName = config().getString("example");
		return ExampleConfigurationFactory.getInstance().createExampleConfiguration(exampleName).net;
	}

}
