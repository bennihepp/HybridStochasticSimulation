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
	private int[][] productionStoichiometry;
	private int[][] consumptionStoichiometry;
	private int[][] stoichiometry;
	private double[] rateParameters;
	private transient volatile List<int[]> choiceIndicesList = null;
	private transient volatile List<List<Integer>> involvedSpeciesList = null;
	private transient volatile List<List<Integer>> involvedReactionsList = null;
	private transient volatile ReactionNetworkGraph graph = null;
	private List<String> speciesLabels;
	private List<String> reactionLabels;

	final transient private Object choiceIndicesListMutex = new Object();
	final transient private Object involvedSpeciesListMutex = new Object();
	final transient private Object involvedReactionsListMutex = new Object();
	final transient private Object graphMutex = new Object();

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
				if (network.getStoichiometry(species, reaction) != 0) {
					subReactions.put(j, reaction);
					j++;
					break;
				}
			}
		}
		int[][] subProductionStoichiometries = new int[subReactions.size()][subnetworkSpecies.size()];
		int[][] subConsumptionStoichiometries = new int[subReactions.size()][subnetworkSpecies.size()];
		DefaultUnaryBinaryReactionNetwork subNetwork = new DefaultUnaryBinaryReactionNetwork(subnetworkSpecies.size(), subReactions.size());
		for (Map.Entry<Integer, Integer> entry : subSpeciesMap.entrySet()) {
			int subSpecies = entry.getKey();
			subNetwork.setSpeciesLabel(subSpecies, network.getSpeciesLabel(entry.getValue()));
		}
		for (Map.Entry<Integer, Integer> entry : subReactions.entrySet()) {
			int subReaction = entry.getKey();
			double rateParameter = network.getRateParameter(entry.getValue());
			subNetwork.setRateParameter(subReaction, rateParameter);
			int[] productionStoichiometries = network.getProductionStoichiometries(entry.getValue());
			int[] consumptionStoichiometries = network.getConsumptionStoichiometries(entry.getValue());;
			for (int k=0; k < subNetwork.getNumberOfSpecies(); k++) {
				subProductionStoichiometries[subReaction][k] = productionStoichiometries[subSpeciesMap.get(k)];
				subConsumptionStoichiometries[subReaction][k] = consumptionStoichiometries[subSpeciesMap.get(k)];
			}
			subNetwork.setReactionLabel(subReaction, network.getReactionLabel(entry.getValue()));
		}
		subNetwork.setStoichiometries(subProductionStoichiometries, subConsumptionStoichiometries);
		subNetwork.graph = SubReactionNetworkGraph.createSubReactionNetworkGraph(network, subVerticesMap, subEdgesMap);
		return subNetwork;
	}

	public DefaultUnaryBinaryReactionNetwork(int numOfSpecies, int numOfReactions) {
		checkArgument(numOfSpecies > 0, "Expected numOfSpecies > 0");
		checkArgument(numOfReactions >= 0, "Expected numOfReactions > 0");
		this.numOfSpecies = numOfSpecies;
		this.numOfReactions = numOfReactions;
		initialConditions = new double[numOfSpecies];
		productionStoichiometry = new int[numOfReactions][numOfSpecies];
		consumptionStoichiometry = new int[numOfReactions][numOfSpecies];
		stoichiometry = new int[numOfReactions][numOfSpecies];
		rateParameters = new double[numOfReactions];
		initLabels();
	}

	protected DefaultUnaryBinaryReactionNetwork(UnaryBinaryReactionNetwork network) {
		this(network.getNumberOfSpecies(), network.getNumberOfReactions());
		setSpeciesLabels(network.getSpeciesLabels());
		setReactionLabels(network.getReactionLabels());
		setRateParameters(network.getRateParameters());
		setStoichiometries(network.getProductionStoichiometries(), network.getConsumptionStoichiometries());
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

	public void setStoichiometry(int species, int reaction, int production, int consumption) {
		checkElementIndex(species, getNumberOfSpecies(), "Expected 0<=species<getNumberOfSpecies()");
		checkElementIndex(reaction, getNumberOfReactions(), "Expected 0<=species<getNumberOfSpecies()");
		checkArgument(production >= 0, "Expected production >= 0");
		checkArgument(consumption >= 0, "Expected consumption >= 0");
		productionStoichiometry[reaction][species] = production;
		consumptionStoichiometry[reaction][species] = consumption;
		stoichiometry[reaction][species] = production - consumption;
		invalidate();
	}

	public void setProductionStoichiometry(int species, int reaction, int production) {
		checkElementIndex(species, getNumberOfSpecies(), "Expected 0<=species<getNumberOfSpecies()");
		checkElementIndex(reaction, getNumberOfReactions(), "Expected 0<=species<getNumberOfSpecies()");
		checkArgument(production >= 0, "Expected production >= 0");
		productionStoichiometry[reaction][species] = production;
		int consumption = consumptionStoichiometry[reaction][species];
		stoichiometry[reaction][species] = production - consumption;
		invalidate();
	}

	public void setConsumptionStoichiometry(int species, int reaction, int consumption) {
		checkElementIndex(species, getNumberOfSpecies(), "Expected 0<=species<getNumberOfSpecies()");
		checkElementIndex(reaction, getNumberOfReactions(), "Expected 0<=species<getNumberOfSpecies()");
		checkArgument(consumption >= 0, "Expected consumption >= 0");
		int production = productionStoichiometry[reaction][species];
		consumptionStoichiometry[reaction][species] = consumption;
		stoichiometry[reaction][species] = production - consumption;
		invalidate();
	}

	public void setStoichiometries(int[][] productionStoich, int[][] consumptionStoich) {
		checkArgument(productionStoich.length == getNumberOfReactions());
		checkArgument(consumptionStoich.length == getNumberOfReactions());
		if (getNumberOfReactions() == 0)
			return;
		checkArgument(productionStoich[0].length == getNumberOfSpecies());
		checkArgument(consumptionStoich[0].length == getNumberOfSpecies());
		for (int r=0; r < getNumberOfReactions(); r++) {
			System.arraycopy(productionStoich[r], 0, productionStoichiometry[r], 0, numOfSpecies);
			System.arraycopy(consumptionStoich[r], 0, consumptionStoichiometry[r], 0, numOfSpecies);
			for (int s=0; s < getNumberOfSpecies(); s++)
				stoichiometry[r][s] = productionStoich[r][s] - consumptionStoich[r][s];
		}
		invalidate();
	}

	@Override
	final public int getProductionStoichiometry(int species, int reaction) {
		checkElementIndex(species, getNumberOfSpecies());
		checkElementIndex(reaction, getNumberOfReactions());
		return productionStoichiometry[reaction][species];
	}

	@Override
	final public int getConsumptionStoichiometry(int species, int reaction) {
		checkElementIndex(species, getNumberOfSpecies());
		checkElementIndex(reaction, getNumberOfReactions());
		return consumptionStoichiometry[reaction][species];
	}

	@Override
	final public int getStoichiometry(int species, int reaction) {
		checkElementIndex(species, getNumberOfSpecies());
		checkElementIndex(reaction, getNumberOfReactions());
		return stoichiometry[reaction][species];
	}

	@Override
	public int[] getProductionStoichiometries(int reaction) {
		checkElementIndex(reaction, getNumberOfReactions());
		return productionStoichiometry[reaction].clone();
	}

	@Override
	public int[] getConsumptionStoichiometries(int reaction) {
		checkElementIndex(reaction, getNumberOfReactions());
		return consumptionStoichiometry[reaction].clone();
	}

	@Override
	public int[] getStoichiometries(int reaction) {
		checkElementIndex(reaction, getNumberOfReactions());
		return stoichiometry[reaction].clone();
	}

	@Override
	public int[][] getProductionStoichiometries() {
		return productionStoichiometry.clone();
	}

	@Override
	public int[][] getConsumptionStoichiometries() {
		return consumptionStoichiometry.clone();
	}

	@Override
	public int[][] getStoichiometries() {
		return stoichiometry.clone();
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
				switch (consumptionStoichiometry[r][s]) {
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
				if (getProductionStoichiometry(s, r) != 0 || getConsumptionStoichiometry(s, r) != 0)
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
				if (getProductionStoichiometry(s, r) != 0 || getConsumptionStoichiometry(s, r) != 0)
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
