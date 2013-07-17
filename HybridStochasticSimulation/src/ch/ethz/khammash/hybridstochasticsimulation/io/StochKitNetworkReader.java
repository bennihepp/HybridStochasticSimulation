package ch.ethz.khammash.hybridstochasticsimulation.io;

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

import ch.ethz.khammash.hybridstochasticsimulation.examples.SimulationConfiguration;
import ch.ethz.khammash.hybridstochasticsimulation.networks.DefaultUnaryBinaryReactionNetwork;
import ch.ethz.khammash.hybridstochasticsimulation.networks.UnaryBinaryReactionNetwork;

public class StochKitNetworkReader {

	private StochKitNetworkReader() {}

	public static UnaryBinaryReactionNetwork readUnaryBinaryNetworkFromFile(File inputFile) throws ParserConfigurationException, SAXException, IOException {
		return readSimulationConfiguration(inputFile).net;
	}

	public static SimulationConfiguration readSimulationConfiguration(File inputFile) throws ParserConfigurationException, SAXException, IOException {
		DocumentBuilderFactory factory = DocumentBuilderFactory.newInstance();
	    DocumentBuilder builder = factory.newDocumentBuilder();
	    Document document = builder.parse(inputFile);

	    NodeList nodeList = document.getElementsByTagName("NumberOfReactions");
	    int numOfReactions = Integer.parseInt(nodeList.item(0).getTextContent());
	    nodeList = document.getElementsByTagName("NumberOfSpecies");
	    int numOfSpecies = Integer.parseInt(nodeList.item(0).getTextContent());

	    DefaultUnaryBinaryReactionNetwork net = new DefaultUnaryBinaryReactionNetwork(numOfSpecies, numOfReactions);

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
	    	parseReactions(net, speciesMap, parameterMap, (Element)nodeList.item(i));

	    double[] x0 = new double[net.getNumberOfSpecies()];
	    String[] labels = new String[net.getNumberOfSpecies()];
	    for (Map.Entry<String, Integer> entry : speciesMap.entrySet()) {
	    	int species = entry.getValue();
	    	double value = initialPopulationMap.get(entry.getKey());
	    	x0[species] = value;
	    	labels[species] = entry.getKey();
	    }

	    SimulationConfiguration nss = new SimulationConfiguration();
	    nss.net = net;
	    nss.reset();
	    nss.x0 = x0;
	    nss.speciesNames = labels;
	    return nss;
	}

	private static void parseReactions(DefaultUnaryBinaryReactionNetwork net,
			Map<String, Integer> speciesMap, Map<String, Double> parameterMap, Element element) {
		NodeList nodeList = element.getElementsByTagName("Reaction");
		for (int reaction=0; reaction < nodeList.getLength(); reaction++) {
			Element childElement = (Element)nodeList.item(reaction);
//			String id = childElement.getElementsByTagName("Id").item(0).getTextContent();
			String type = childElement.getElementsByTagName("Type").item(0).getTextContent();
			if (!type.equals("mass-action"))
				// TODO: Custom exception type
				throw new RuntimeException("Unsupported reaction type: " + type);
			String rateString = childElement.getElementsByTagName("Rate").item(0).getTextContent();
			double rate = parameterMap.get(rateString);
			net.setRateParameter(reaction, rate);
			Map<String,Integer> reactantMap = readReactants(childElement);
			for (Entry<String, Integer> entry : reactantMap.entrySet()) {
				int species = speciesMap.get(entry.getKey());
				int consumption = entry.getValue();
				net.setConsumptionStochiometry(species, reaction, consumption);
			}
			Map<String,Integer> productMap = readProducts(childElement);
			for (Entry<String, Integer> entry : productMap.entrySet()) {
				int species = speciesMap.get(entry.getKey());
				int production = entry.getValue();
				net.setProductionStochiometry(species, reaction, production);
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
