package ch.ethz.bhepp.hybridstochasticsimulation.io;

import java.io.File;
import java.io.IOException;
import java.util.HashMap;
import java.util.Map;
import java.util.Map.Entry;

import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;
import javax.xml.parsers.ParserConfigurationException;

import org.apache.commons.math3.util.FastMath;
import org.w3c.dom.Document;
import org.w3c.dom.Element;
import org.w3c.dom.NodeList;
import org.xml.sax.SAXException;

import ch.ethz.bhepp.hybridstochasticsimulation.examples.SimulationConfiguration;
import ch.ethz.bhepp.hybridstochasticsimulation.networks.DefaultUnaryBinaryReactionNetwork;
import ch.ethz.bhepp.hybridstochasticsimulation.networks.MassActionReactionNetwork;

public class StochKitNetworkReader {

	public static class FileFormatException extends Exception {
		private static final long serialVersionUID = -8291451734049197094L;

		public FileFormatException(String message) {
			super(message);
		}
	}

	public static class NoSuchParameterException extends FileFormatException {

		private static final long serialVersionUID = 4373584458629447987L;

		public NoSuchParameterException(String message) {
			super(message);
		}
	}

	private StochKitNetworkReader() {}

	public static MassActionReactionNetwork readUnaryBinaryNetworkFromFile(File inputFile) throws ParserConfigurationException, SAXException, IOException, FileFormatException {
		return readSimulationConfiguration(inputFile).net;
	}

	public static SimulationConfiguration readSimulationConfiguration(File inputFile) throws ParserConfigurationException, SAXException, IOException, FileFormatException {
		DocumentBuilderFactory factory = DocumentBuilderFactory.newInstance();
	    DocumentBuilder builder = factory.newDocumentBuilder();
	    Document document = builder.parse(inputFile);

	    int numOfReactions;
	    int numOfSpecies;
	    NodeList nodeList = document.getElementsByTagName("NumberOfReactions");
	    if (nodeList.getLength() > 0)
	    	numOfReactions = Integer.parseInt(nodeList.item(0).getTextContent());
	    else
	    	numOfReactions = document.getElementsByTagName("Reaction").getLength();
	    nodeList = document.getElementsByTagName("NumberOfSpecies");
	    if (nodeList.getLength() > 0)
		    numOfSpecies = Integer.parseInt(nodeList.item(0).getTextContent());
	    else
	    	numOfSpecies = document.getElementsByTagName("Species").getLength();

	    DefaultUnaryBinaryReactionNetwork network = new DefaultUnaryBinaryReactionNetwork(numOfSpecies, numOfReactions);

	    HashMap<String, Double> parameterMap = new HashMap<>();
	    nodeList = document.getElementsByTagName("ParametersList");
	    for (int i=0; i < nodeList.getLength(); i++)
	    	parameterMap.putAll(readParameters((Element)nodeList.item(i)));

	    HashMap<String, Double> initialPopulationMap = new HashMap<>();
	    HashMap<String, Integer> speciesMap = new HashMap<>();
	    nodeList = document.getElementsByTagName("SpeciesList");
	    for (int i=0; i < nodeList.getLength(); i++)
	    	initialPopulationMap.putAll(readSpecies((Element)nodeList.item(i)));
	    for (String key : initialPopulationMap.keySet())
	    	speciesMap.put(key, speciesMap.size());

	    nodeList = document.getElementsByTagName("ReactionsList");
	    for (int i=0; i < nodeList.getLength(); i++)
	    	parseReactions(network, speciesMap, parameterMap, (Element)nodeList.item(i));

	    double[] x0 = new double[network.getNumberOfSpecies()];
	    String[] labels = new String[network.getNumberOfSpecies()];
	    for (Map.Entry<String, Integer> entry : speciesMap.entrySet()) {
	    	int species = entry.getValue();
	    	double value = initialPopulationMap.get(entry.getKey());
	    	x0[species] = value;
	    	network.setSpeciesLabel(species, entry.getKey());
	    	labels[species] = entry.getKey();
	    }

	    network.setInitialConditions(x0);

	    SimulationConfiguration nss = new SimulationConfiguration();
	    nss.net = network;
	    nss.reset();
	    nss.x0 = x0;
	    nss.speciesNames = labels;
	    return nss;
	}

	private static void parseReactions(DefaultUnaryBinaryReactionNetwork network,
			Map<String, Integer> speciesMap, Map<String, Double> parameterMap, Element element) throws FileFormatException {
		NodeList nodeList = element.getElementsByTagName("Reaction");
		for (int reaction=0; reaction < nodeList.getLength(); reaction++) {
			Element childElement = (Element)nodeList.item(reaction);
			String id = childElement.getElementsByTagName("Id").item(0).getTextContent();
	    	network.setReactionLabel(reaction, id);
			String type = childElement.getElementsByTagName("Type").item(0).getTextContent();
			if (!type.equals("mass-action"))
				throw new FileFormatException("Unsupported reaction type: " + type);
			String rateString = childElement.getElementsByTagName("Rate").item(0).getTextContent();
			if (!parameterMap.containsKey(rateString))
				throw new NoSuchParameterException(rateString);
			double rate = parameterMap.get(rateString);
			network.setRateParameter(reaction, rate);
			Map<String,Integer> reactantMap = readReactants(childElement);
			for (Entry<String, Integer> entry : reactantMap.entrySet()) {
				int species = speciesMap.get(entry.getKey());
				int consumption = entry.getValue();
				network.setConsumptionStoichiometry(species, reaction, consumption);
			}
			Map<String,Integer> productMap = readProducts(childElement);
			for (Entry<String, Integer> entry : productMap.entrySet()) {
				int species = speciesMap.get(entry.getKey());
				int production = entry.getValue();
				network.setProductionStoichiometry(species, reaction, production);
			}
		}
	}

	private static Map<String, Integer> readProducts(Element element) {
		HashMap<String, Integer> productMap = new HashMap<String, Integer>();
		NodeList nodeList = element.getElementsByTagName("Products");
		if (nodeList.getLength() <= 0)
			return productMap;
		Element productsElement = (Element)nodeList.item(0);
		nodeList = productsElement.getElementsByTagName("SpeciesReference");
		for (int i=0; i < nodeList.getLength(); i++) {
			Element childElement = (Element)nodeList.item(i);
			String id = childElement.getAttribute("id");
			int stochiometry = (int)FastMath.round(Double.parseDouble(childElement.getAttribute("stoichiometry")));
			productMap.put(id, stochiometry);
		}
		return productMap;
	}

	private static Map<String, Integer> readReactants(Element element) {
		HashMap<String, Integer> reactantMap = new HashMap<String, Integer>();
		NodeList nodeList = element.getElementsByTagName("Reactants");
		if (nodeList.getLength() <= 0)
			return reactantMap;
		Element reactantsElement = (Element)nodeList.item(0);
		nodeList = reactantsElement.getElementsByTagName("SpeciesReference");
		for (int i=0; i < nodeList.getLength(); i++) {
			Element childElement = (Element)nodeList.item(i);
			String id = childElement.getAttribute("id");
			int stochiometry = (int)FastMath.round(Double.parseDouble(childElement.getAttribute("stoichiometry")));
			reactantMap.put(id, stochiometry);
		}
		return reactantMap;
	}

	private static HashMap<String, Double> readSpecies(Element element) {
		HashMap<String, Double> speciesMap = new HashMap<String, Double>();
		NodeList nodeList = element.getElementsByTagName("Species");
		for (int i=0; i < nodeList.getLength(); i++) {
			Element childElement = (Element)nodeList.item(i);
			String id = childElement.getElementsByTagName("Id").item(0).getTextContent();
			Double value = Double.parseDouble(childElement.getElementsByTagName("InitialPopulation").item(0).getTextContent());
			speciesMap.put(id, value);
		}
		return speciesMap;
	}

	private static HashMap<String, Double> readParameters(Element element) {
		HashMap<String, Double> parameterMap = new HashMap<String, Double>();
		NodeList nodeList = element.getElementsByTagName("Parameter");
		for (int i=0; i < nodeList.getLength(); i++) {
			Element childElement = (Element)nodeList.item(i);
			String id = childElement.getElementsByTagName("Id").item(0).getTextContent();
			Double value = Double.parseDouble(childElement.getElementsByTagName("Expression").item(0).getTextContent());
			parameterMap.put(id, value);
		}
		return parameterMap;
	}
}
