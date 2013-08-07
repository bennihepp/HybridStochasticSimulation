package ch.ethz.khammash.hybridstochasticsimulation.networks;

import static com.google.common.base.Preconditions.checkArgument;
import static com.google.common.base.Preconditions.checkElementIndex;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;

import ch.ethz.khammash.hybridstochasticsimulation.graphs.DefaultReactionNetworkGraph;
import ch.ethz.khammash.hybridstochasticsimulation.graphs.ReactionEdge;
import ch.ethz.khammash.hybridstochasticsimulation.graphs.ReactionNetworkGraph;
import ch.ethz.khammash.hybridstochasticsimulation.graphs.SpeciesVertex;
import ch.ethz.khammash.hybridstochasticsimulation.graphs.SubReactionNetworkGraph;


/*
 * rateParameters k_i are of the form so that the propensity of a reaction is
 * 	k_i for constitutive reactions
 *  k_i * x1 for unary reactions
 *  k_i * x1 * x2 for binary reactions of different species
 *  1/2 * k_i * x1 * (x1 - 1) for binary reactions of the same species
 */
public class DefaultUnaryBinaryReactionNetwork implements UnaryBinaryReactionNetwork {

	private int numOfSpecies;
	private int numOfReactions;
	private double[] initialConditions;
	private int[][] productionStochiometry;
	private int[][] consumptionStochiometry;
	private int[][] stochiometry;
	private double[] rateParameters;
	private volatile List<int[]> choiceIndicesList;
	private volatile List<List<Integer>> involvedSpeciesList;
	private volatile List<List<Integer>> involvedReactionsList;
	private volatile ReactionNetworkGraph graph;
	private List<String> speciesLabels;
	private List<String> reactionLabels;

	final private Object choiceIndicesListMutex = new Object();
	final private Object involvedSpeciesListMutex = new Object();
	final private Object involvedReactionsListMutex = new Object();
	final private Object graphMutex = new Object();

	public static DefaultUnaryBinaryReactionNetwork createFromNetwork(UnaryBinaryReactionNetwork network) {
		return new DefaultUnaryBinaryReactionNetwork(network);
	}

	public static DefaultUnaryBinaryReactionNetwork createSubnetwork(UnaryBinaryReactionNetwork network, Set<Integer> subnetworkSpecies) {
		Map<Integer, SpeciesVertex> subVerticesMap = new HashMap<>();
		Map<Integer, List<ReactionEdge>> subEdgesMap = new HashMap<>();
		Map<Integer, Integer> subSpeciesMap = new HashMap<>(subnetworkSpecies.size());
		int i = 0;
		for (int species : subnetworkSpecies) {
			subSpeciesMap.put(i, species);
			subVerticesMap.put(i, network.getGraph().getSpeciesVertex(species));
			i++;
		}
		Map<Integer, Integer> subReactions = new HashMap<>(network.getNumberOfReactions());
		int j = 0;
		for (int reaction=0; reaction < network.getNumberOfReactions(); reaction++) {
			boolean interactsWithOutsideSpecies = false;
			List<Integer> involvedSpecies = network.getInvolvedSpecies(reaction);
			for (int species : involvedSpecies)
				if (!subnetworkSpecies.contains(species)) {
					interactsWithOutsideSpecies = true;
					break;
				}
			if (interactsWithOutsideSpecies)
				continue;
			List<ReactionEdge> edgeList = new LinkedList<>();
			for (ReactionEdge edge : network.getGraph().getReactionEdges(reaction)) {
				if (subnetworkSpecies.contains(edge.getSource().getSpecies()) &&
						subnetworkSpecies.contains(edge.getTarget().getSpecies()))
						edgeList.add(edge);
			}
			subEdgesMap.put(j, edgeList);
			for (int species : subnetworkSpecies) {
				if (network.getStochiometry(species, reaction) != 0) {
					subReactions.put(j, reaction);
					j++;
					break;
				}
			}
		}
		int[][] subProductionStochiometries = new int[subReactions.size()][subnetworkSpecies.size()];
		int[][] subConsumptionStochiometries = new int[subReactions.size()][subnetworkSpecies.size()];
		DefaultUnaryBinaryReactionNetwork subNetwork = new DefaultUnaryBinaryReactionNetwork(subnetworkSpecies.size(), subReactions.size());
		for (Map.Entry<Integer, Integer> entry : subSpeciesMap.entrySet()) {
			int subSpecies = entry.getKey();
			subNetwork.setSpeciesLabel(subSpecies, network.getSpeciesLabel(entry.getValue()));
		}
		for (Map.Entry<Integer, Integer> entry : subReactions.entrySet()) {
			int subReaction = entry.getKey();
			double rateParameter = network.getRateParameter(entry.getValue());
			subNetwork.setRateParameter(subReaction, rateParameter);
			int[] productionStochiometries = network.getProductionStochiometries(entry.getValue());
			int[] consumptionStochiometries = network.getConsumptionStochiometries(entry.getValue());;
			for (int k=0; k < subNetwork.getNumberOfSpecies(); k++) {
				subProductionStochiometries[subReaction][k] = productionStochiometries[subSpeciesMap.get(k)];
				subConsumptionStochiometries[subReaction][k] = consumptionStochiometries[subSpeciesMap.get(k)];
			}
			subNetwork.setReactionLabel(subReaction, network.getReactionLabel(entry.getValue()));
		}
		subNetwork.setStochiometries(subProductionStochiometries, subConsumptionStochiometries);
		subNetwork.graph = SubReactionNetworkGraph.createSubReactionNetworkGraph(network, subVerticesMap, subEdgesMap);
		return subNetwork;
	}

	public DefaultUnaryBinaryReactionNetwork(int numOfSpecies, int numOfReactions) {
		checkArgument(numOfSpecies > 0, "Expected numOfSpecies > 0");
		checkArgument(numOfReactions >= 0, "Expected numOfReactions > 0");
		this.numOfSpecies = numOfSpecies;
		this.numOfReactions = numOfReactions;
		initialConditions = new double[numOfSpecies];
		productionStochiometry = new int[numOfReactions][numOfSpecies];
		consumptionStochiometry = new int[numOfReactions][numOfSpecies];
		stochiometry = new int[numOfReactions][numOfSpecies];
		rateParameters = new double[numOfReactions];
		initLabels();
	}

	protected DefaultUnaryBinaryReactionNetwork(UnaryBinaryReactionNetwork network) {
		this(network.getNumberOfSpecies(), network.getNumberOfReactions());
		setSpeciesLabels(network.getSpeciesLabels());
		setReactionLabels(network.getReactionLabels());
		setRateParameters(network.getRateParameters());
		setStochiometries(network.getProductionStochiometries(), network.getConsumptionStochiometries());
		setInitialConditions(network.getInitialConditions());
	}

	final private void initLabels() {
		speciesLabels = new ArrayList<>(numOfSpecies);
		for (int s=0; s < numOfSpecies; s++)
			speciesLabels.add("S" + s);
		reactionLabels = new ArrayList<>(numOfReactions);
		for (int r=0; r < numOfReactions; r++)
			reactionLabels.add("R" + r);
	}

	@Override
	final public int getNumberOfSpecies() {
		return numOfSpecies;
	}

	@Override
	final public int getNumberOfReactions() {
		return numOfReactions;
	}

	final protected void invalidate() {
		choiceIndicesList = null;
		involvedSpeciesList = null;
		involvedReactionsList = null;
	}

	public void setStochiometry(int species, int reaction, int production, int consumption) {
		checkElementIndex(species, getNumberOfSpecies(), "Expected 0<=species<getNumberOfSpecies()");
		checkElementIndex(reaction, getNumberOfReactions(), "Expected 0<=species<getNumberOfSpecies()");
		checkArgument(production >= 0, "Expected production >= 0");
		checkArgument(consumption >= 0, "Expected consumption >= 0");
		productionStochiometry[reaction][species] = production;
		consumptionStochiometry[reaction][species] = consumption;
		stochiometry[reaction][species] = production - consumption;
		invalidate();
	}

	public void setProductionStochiometry(int species, int reaction, int production) {
		checkElementIndex(species, getNumberOfSpecies(), "Expected 0<=species<getNumberOfSpecies()");
		checkElementIndex(reaction, getNumberOfReactions(), "Expected 0<=species<getNumberOfSpecies()");
		checkArgument(production >= 0, "Expected production >= 0");
		productionStochiometry[reaction][species] = production;
		int consumption = consumptionStochiometry[reaction][species];
		stochiometry[reaction][species] = production - consumption;
		invalidate();
	}

	public void setConsumptionStochiometry(int species, int reaction, int consumption) {
		checkElementIndex(species, getNumberOfSpecies(), "Expected 0<=species<getNumberOfSpecies()");
		checkElementIndex(reaction, getNumberOfReactions(), "Expected 0<=species<getNumberOfSpecies()");
		checkArgument(consumption >= 0, "Expected consumption >= 0");
		int production = productionStochiometry[reaction][species];
		consumptionStochiometry[reaction][species] = consumption;
		stochiometry[reaction][species] = production - consumption;
		invalidate();
	}

	public void setStochiometries(int[][] productionStoch, int[][] consumptionStoch) {
		checkArgument(productionStoch.length == getNumberOfReactions());
		checkArgument(consumptionStoch.length == getNumberOfReactions());
		if (getNumberOfReactions() == 0)
			return;
		checkArgument(productionStoch[0].length == getNumberOfSpecies());
		checkArgument(consumptionStoch[0].length == getNumberOfSpecies());
		for (int r=0; r < getNumberOfReactions(); r++) {
			System.arraycopy(productionStoch[r], 0, productionStochiometry[r], 0, numOfSpecies);
			System.arraycopy(consumptionStoch[r], 0, consumptionStochiometry[r], 0, numOfSpecies);
			for (int s=0; s < getNumberOfSpecies(); s++)
				stochiometry[r][s] = productionStoch[r][s] - consumptionStoch[r][s];
		}
		invalidate();
	}

	@Override
	final public int getProductionStochiometry(int species, int reaction) {
		checkElementIndex(species, getNumberOfSpecies());
		checkElementIndex(reaction, getNumberOfReactions());
		return productionStochiometry[reaction][species];
	}

	@Override
	final public int getConsumptionStochiometry(int species, int reaction) {
		checkElementIndex(species, getNumberOfSpecies());
		checkElementIndex(reaction, getNumberOfReactions());
		return consumptionStochiometry[reaction][species];
	}

	@Override
	final public int getStochiometry(int species, int reaction) {
		checkElementIndex(species, getNumberOfSpecies());
		checkElementIndex(reaction, getNumberOfReactions());
		return stochiometry[reaction][species];
	}

	@Override
	public int[] getProductionStochiometries(int reaction) {
		checkElementIndex(reaction, getNumberOfReactions());
		return productionStochiometry[reaction].clone();
	}

	@Override
	public int[] getConsumptionStochiometries(int reaction) {
		checkElementIndex(reaction, getNumberOfReactions());
		return consumptionStochiometry[reaction].clone();
	}

	@Override
	public int[] getStochiometries(int reaction) {
		checkElementIndex(reaction, getNumberOfReactions());
		return stochiometry[reaction].clone();
	}

	@Override
	public int[][] getProductionStochiometries() {
		return productionStochiometry.clone();
	}

	@Override
	public int[][] getConsumptionStochiometries() {
		return consumptionStochiometry.clone();
	}

	@Override
	public int[][] getStochiometries() {
		return stochiometry.clone();
	}

	public void setRateParameter(int reaction, double parameter) {
		checkElementIndex(reaction, getNumberOfReactions());
		this.rateParameters[reaction] = parameter;
	}

	public void setRateParameters(double[] rateParameters) {
		System.arraycopy(rateParameters, 0, this.rateParameters, 0, getNumberOfReactions());
		for (int r=0; r < numOfReactions; r++)
			this.rateParameters[r] = rateParameters[r];
	}

	@Override
	final public double getRateParameter(int reaction) {
		checkElementIndex(reaction, getNumberOfReactions());
		return rateParameters[reaction];
	}

	@Override
	public double[] getRateParameters() {
		return rateParameters.clone();
	}

	// Returns an array of int, where each entry is a consumed species of the reaction
	@Override
	public int[] getChoiceIndices(int reaction) {
		checkElementIndex(reaction, getNumberOfReactions());
		return getModifieableChoiceIndicesList().get(reaction);
	}

	@Override
	public List<int[]> getChoiceIndicesList() {
		return Collections.unmodifiableList(getModifieableChoiceIndicesList());
	}

	private List<int[]> getModifieableChoiceIndicesList() {
		if (choiceIndicesList == null) {
			synchronized (choiceIndicesListMutex) {
				if (choiceIndicesList == null) {
					choiceIndicesList = computeChoiceIndices();
				}
			}
		}
		return choiceIndicesList;
	}

	private List<int[]> computeChoiceIndices() {
		List<int[]> choiceIndicesList = new ArrayList<>();
		for (int r=0; r < getNumberOfReactions(); r++) {
			int s1 = -1;
			int s2 = -1;
			for (int s=0; s < getNumberOfSpecies(); s++)
				switch (consumptionStochiometry[r][s]) {
				case 0:
					break;
				case 1:
					if (s1 == -1)
						s1 = s;
					else if (s2 == -1)
						s2 = s;
					else
						throw new RuntimeException("Only constitutive, unary and binary reactions are allowed");
					break;
				case 2:
					if (s1 == -1) {
						s1 = s;
						s2 = s;
					} else
						throw new RuntimeException("Only constitutive, unary and binary reactions are allowed");
					break;
				default:
					throw new RuntimeException("Only constitutive, unary and binary reactions are allowed");
				}
			int[] choiceIndices = null;
			if (s1 != -1 && s2 != -1) {
				choiceIndices = new int[2];
				choiceIndices[0] = s1;
				choiceIndices[1] = s2;
			} else if (s1 != -1) {
				choiceIndices = new int[1];
				choiceIndices[0] = s1;
			} else if (s2 != -1) {
				choiceIndices = new int[1];
				choiceIndices[0] = s2;
			} else {
				choiceIndices = new int[0];
			}
			choiceIndicesList.add(choiceIndices);
		}
		return choiceIndicesList;
	}

	@Override
	public List<Integer> getInvolvedSpecies(int reaction) {
		return getInvolvedSpeciesList().get(reaction);
	}

	private List<List<Integer>> getInvolvedSpeciesList() {
		if (involvedSpeciesList == null) {
			synchronized (involvedSpeciesListMutex) {
				if (involvedSpeciesList == null) {
					involvedSpeciesList = computeInvolvedSpeciesList();
				}
			}
		}
		return involvedSpeciesList;
	}

	private List<List<Integer>> computeInvolvedSpeciesList() {
		List<List<Integer>> involvedSpecies = new ArrayList<>(getNumberOfReactions());
		for (int r=0; r < getNumberOfReactions(); r++) {
			ArrayList<Integer> is = new ArrayList<Integer>();
			for (int s=0; s < getNumberOfSpecies(); s++)
				if (getProductionStochiometry(s, r) != 0 || getConsumptionStochiometry(s, r) != 0)
					is.add(s);
			involvedSpecies.add(Collections.unmodifiableList(is));
		}
		return involvedSpecies;
	}

	@Override
	public List<Integer> getInvolvedReactions(int species) {
		return getInvolvedReactionsList().get(species);
	}

	private List<List<Integer>> getInvolvedReactionsList() {
		if (involvedReactionsList == null) {
			synchronized (involvedReactionsListMutex) {
				if (involvedReactionsList == null) {
					involvedReactionsList = computeInvolvedReactionsList();
				}
			}
		}
		return involvedReactionsList;
	}

	private List<List<Integer>> computeInvolvedReactionsList() {
		List<List<Integer>> involvedReactionsList = new ArrayList<>(getNumberOfSpecies());
		for (int s=0; s < getNumberOfSpecies(); s++) {
			ArrayList<Integer> ir = new ArrayList<Integer>();
			for (int r=0; r < getNumberOfReactions(); r++)
				if (getProductionStochiometry(s, r) != 0 || getConsumptionStochiometry(s, r) != 0)
					ir.add(r);
			involvedReactionsList.add(Collections.unmodifiableList(ir));
		}
		return involvedReactionsList;
	}

	@Override
	public ReactionNetworkGraph getGraph() {
		if (graph == null) {
			synchronized (graphMutex) {
				if (graph == null) {
					graph = new DefaultReactionNetworkGraph(this);
				}
			}
		}
		return graph;
	}

	@Override
	public String getSpeciesLabel(int species) {
		return speciesLabels.get(species);
	}

	@Override
	public String getReactionLabel(int reaction) {
		return reactionLabels.get(reaction);
	}

	@Override
	public List<String> getSpeciesLabels() {
		return Collections.unmodifiableList(speciesLabels);
	}

	@Override
	public List<String> getReactionLabels() {
		return Collections.unmodifiableList(reactionLabels);
	}

	public void setSpeciesLabel(int species, String label) {
		speciesLabels.set(species, label);
	}

	public void setReactionLabel(int reaction, String label) {
		reactionLabels.set(reaction, label);
	}

	public void setSpeciesLabels(List<String> speciesLabels) {
		this.speciesLabels = new ArrayList<>(speciesLabels);
	}

	public void setReactionLabels(List<String> reactionLabels) {
		this.reactionLabels = new ArrayList<>(reactionLabels);
	}

	@Override
	public double[] getInitialConditions() {
		return initialConditions.clone();
	}

	public void setInitialConditions(double[] initialConditions) {
		checkArgument(initialConditions.length == getNumberOfSpecies(), "Expected initialConditions.length == getNumberOfSpecies()");
		System.arraycopy(initialConditions, 0, this.initialConditions, 0, getNumberOfSpecies());
		this.initialConditions = initialConditions;
	}

}
