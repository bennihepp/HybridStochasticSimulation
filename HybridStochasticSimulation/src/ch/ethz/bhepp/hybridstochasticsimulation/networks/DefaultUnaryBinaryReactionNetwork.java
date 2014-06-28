package ch.ethz.bhepp.hybridstochasticsimulation.networks;

import static com.google.common.base.Preconditions.checkArgument;
import static com.google.common.base.Preconditions.checkElementIndex;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.apache.commons.lang3.ArrayUtils;

import ch.ethz.bhepp.hybridstochasticsimulation.graphs.DefaultReactionNetworkGraph;
import ch.ethz.bhepp.hybridstochasticsimulation.graphs.ReactionEdge;
import ch.ethz.bhepp.hybridstochasticsimulation.graphs.ReactionNetworkGraph;
import ch.ethz.bhepp.hybridstochasticsimulation.graphs.SpeciesVertex;
import ch.ethz.bhepp.hybridstochasticsimulation.graphs.SubReactionNetworkGraph;

import com.google.common.base.Supplier;
import com.google.common.base.Suppliers;
import com.google.common.collect.ImmutableSet;


/*
 * rateParameters k_i are of the form so that the propensity of a reaction is
 * 	k_i for constitutive reactions
 *  k_i * x1 for unary reactions
 *  k_i * x1 * x2 for binary reactions of different species
 *  1/2 * k_i * x1 * (x1 - 1) for binary reactions of the same species
 */
public class DefaultUnaryBinaryReactionNetwork implements MassActionReactionNetwork {

	private int numOfSpecies;
	private int numOfReactions;
	private double[] initialConditions;
	private int[][] productionStoichiometry;
	private int[][] consumptionStoichiometry;
	private int[][] stoichiometry;
	private double[] rateParameters;
	private List<String> speciesLabels;
	private List<String> reactionLabels;

	private Supplier<ReactionNetworkGraph> graphSupplier;
	private Supplier<List<int[]>> reactantIndicesArraySupplier;
	private Supplier<List<List<Integer>>> reactantIndicesListSupplier;
	private Supplier<List<Set<Integer>>> involvedSpeciesListSupplier;
	private Supplier<List<Set<Integer>>> involvingReactionsListSupplier;
	private Supplier<List<Set<Integer>>> modifiedSpeciesListSupplier;
	private Supplier<List<Set<Integer>>> modifyingReactionsListSupplier;

	public static DefaultUnaryBinaryReactionNetwork createFromNetwork(MassActionReactionNetwork network) {
		return new DefaultUnaryBinaryReactionNetwork(network);
	}

	public static DefaultUnaryBinaryReactionNetwork createSubnetwork(final MassActionReactionNetwork network,
			Set<Integer> subnetworkSpecies, Set<Integer> subnetworkReactions) {
		final Map<Integer, SpeciesVertex> subVerticesMap = new HashMap<>();
		final Map<Integer, List<ReactionEdge>> subEdgesMap = new HashMap<>();
		Map<Integer, Integer> subSpeciesMap = new HashMap<>(subnetworkSpecies.size());
		int i = 0;
		for (int species : subnetworkSpecies) {
			subSpeciesMap.put(i, species);
			subVerticesMap.put(i, network.getGraph().getSpeciesVertex(species));
			i++;
		}
		Map<Integer, Integer> subReactions = new HashMap<>(network.getNumberOfReactions());
		int j = 0;
		for (int reaction : subnetworkReactions) {
			List<ReactionEdge> edgeList = new LinkedList<>();
			for (ReactionEdge edge : network.getGraph().getReactionEdges(reaction)) {
				if (subnetworkSpecies.contains(edge.getSource().getSpecies()) &&
						subnetworkSpecies.contains(edge.getTarget().getSpecies()))
						edgeList.add(edge);
			}
			subEdgesMap.put(j, edgeList);
			subReactions.put(j, reaction);
			j++;
		}
//		for (int reaction=0; reaction < network.getNumberOfReactions(); reaction++) {
//			boolean interactsWithOutsideSpecies = false;
//			Set<Integer> involvedSpecies = network.getInvolvedSpecies(reaction);
//			for (int species : involvedSpecies)
//				if (!subnetworkSpecies.contains(species)) {
//					interactsWithOutsideSpecies = true;
//					break;
//				}
//			if (interactsWithOutsideSpecies)
//				continue;
//			List<ReactionEdge> edgeList = new LinkedList<>();
//			for (ReactionEdge edge : network.getGraph().getReactionEdges(reaction)) {
//				if (subnetworkSpecies.contains(edge.getSource().getSpecies()) &&
//						subnetworkSpecies.contains(edge.getTarget().getSpecies()))
//						edgeList.add(edge);
//			}
//			subEdgesMap.put(j, edgeList);
//			for (int species : subnetworkSpecies) {
//				if (network.getStoichiometry(species, reaction) != 0) {
//					subReactions.put(j, reaction);
//					j++;
//					break;
//				}
//			}
//		}
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
		subNetwork.graphSupplier = Suppliers.memoize(new Supplier<ReactionNetworkGraph>() {

			@Override
			public ReactionNetworkGraph get() {
				return SubReactionNetworkGraph.createSubReactionNetworkGraph(network, subVerticesMap, subEdgesMap);
			}

		});
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
		initSuppliers();
	}

	protected DefaultUnaryBinaryReactionNetwork(MassActionReactionNetwork network) {
		this(network.getNumberOfSpecies(), network.getNumberOfReactions());
		setSpeciesLabels(network.getSpeciesLabels());
		setReactionLabels(network.getReactionLabels());
		setRateParameters(network.getRateParameters());
		setStoichiometries(network.getProductionStoichiometries(), network.getConsumptionStoichiometries());
		setInitialConditions(network.getInitialConditions());
		initSuppliers();
	}

	private void initLabels() {
		speciesLabels = new ArrayList<>(numOfSpecies);
		for (int s=0; s < numOfSpecies; s++)
			speciesLabels.add("S" + s);
		reactionLabels = new ArrayList<>(numOfReactions);
		for (int r=0; r < numOfReactions; r++)
			reactionLabels.add("R" + r);
	}

	private void initSuppliers() {
		graphSupplier = Suppliers.memoize(new Supplier<ReactionNetworkGraph>() {

			@Override
			public ReactionNetworkGraph get() {
				return new DefaultReactionNetworkGraph(DefaultUnaryBinaryReactionNetwork.this);
			}

		});
		reactantIndicesArraySupplier = Suppliers.memoize(new Supplier<List<int[]>>() {

			@Override
			public List<int[]> get() {
				return computeReactantIndicesArrays();
			}

		});
		reactantIndicesListSupplier = Suppliers.memoize(new Supplier<List<List<Integer>>>() {

			@Override
			public List<List<Integer>> get() {
				return computeReactantIndicesLists();
			}

		});
		involvedSpeciesListSupplier = Suppliers.memoize(new Supplier<List<Set<Integer>>>() {

			@Override
			public List<Set<Integer>> get() {
				return computeInvolvedSpeciesList();
			}

		});
		involvingReactionsListSupplier = Suppliers.memoize(new Supplier<List<Set<Integer>>>() {

			@Override
			public List<Set<Integer>> get() {
				return computeInvolvingReactionsList();
			}

		});
		modifiedSpeciesListSupplier = Suppliers.memoize(new Supplier<List<Set<Integer>>>() {

			@Override
			public List<Set<Integer>> get() {
				return computeModifiedSpeciesList();
			}

		});
		modifyingReactionsListSupplier = Suppliers.memoize(new Supplier<List<Set<Integer>>>() {

			@Override
			public List<Set<Integer>> get() {
				return computeModifyingReactionsList();
			}

		});
	}

	@Override
	public int getNumberOfSpecies() {
		return numOfSpecies;
	}

	@Override
	public int getNumberOfReactions() {
		return numOfReactions;
	}

	protected void invalidate() {
		initSuppliers();
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

	@Override
	public int[] getReactantIndices(int r) {
//		checkElementIndex(r, getNumberOfReactions());
		return reactantIndicesArraySupplier.get().get(r);
	}

	@Override
	public List<Integer> getReactantIndicesAsList(int r) {
//		checkElementIndex(r, getNumberOfReactions());
		return reactantIndicesListSupplier.get().get(r);
	}

//	private List<int[]> getReactantIndicesArrayList() {
//		if (reactantIndicesArrayList == null) {
//			synchronized (reactantIndicesMutex) {
//				if (reactantIndicesArrayList == null) {
//					reactantIndicesArrayList = computeReactantIndicesArrays();
//				}
//				reactantIndicesListList = new ArrayList<>(getNumberOfReactions());
//				for (int r=0; r < getNumberOfReactions(); r++) {
//					int[] reactantIndicesArray = reactantIndicesArrayList.get(r);
//					List<Integer> reactantList = Arrays.asList(ArrayUtils.toObject(reactantIndicesArray));
//					reactantIndicesListList.add(reactantList);
//				}
//			}
//		}
//		return reactantIndicesArrayList;
//	}

	private List<int[]> computeReactantIndicesArrays() {
		List<int[]> reactantIndicesArrayList = new ArrayList<>(getNumberOfReactions());
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
			int[] reactantIndices = null;
			if (s1 != -1 && s2 != -1) {
				reactantIndices = new int[2];
				reactantIndices[0] = s1;
				reactantIndices[1] = s2;
			} else if (s1 != -1) {
				reactantIndices = new int[1];
				reactantIndices[0] = s1;
			} else if (s2 != -1) {
				reactantIndices = new int[1];
				reactantIndices[0] = s2;
			} else {
				reactantIndices = new int[0];
			}
			reactantIndicesArrayList.add(reactantIndices);
		}
		return reactantIndicesArrayList;
	}

	private List<List<Integer>> computeReactantIndicesLists() {
		List<List<Integer>> reactantIndicesListList = new ArrayList<>(getNumberOfReactions());
		for (int r=0; r < getNumberOfReactions(); r++) {
			int[] reactantIndicesArray = reactantIndicesArraySupplier.get().get(r);
			List<Integer> reactantList = Arrays.asList(ArrayUtils.toObject(reactantIndicesArray));
			reactantIndicesListList.add(reactantList);
		}
		return reactantIndicesListList;
	}

	@Override
	public Set<Integer> getInvolvedSpecies(int r) {
		return involvedSpeciesListSupplier.get().get(r);
	}

	private List<Set<Integer>> computeInvolvedSpeciesList() {
		List<Set<Integer>> involvedSpeciesList = new ArrayList<>(getNumberOfSpecies());
		for (int r=0; r < getNumberOfReactions(); r++) {
			Set<Integer> is = new HashSet<>();
			for (int s=0; s < getNumberOfSpecies(); s++)
				if (getProductionStoichiometry(s, r) != 0 || getConsumptionStoichiometry(s, r) != 0)
					is.add(s);
			involvedSpeciesList.add(ImmutableSet.copyOf(is));
		}
		return involvedSpeciesList;
	}

	@Override
	public Set<Integer> getInvolvingReactions(int s) {
		return involvingReactionsListSupplier.get().get(s);
	}

	private List<Set<Integer>> computeInvolvingReactionsList() {
		List<Set<Integer>> involvingReactionsList = new ArrayList<>(getNumberOfReactions());
		for (int s=0; s < getNumberOfSpecies(); s++) {
			Set<Integer> ir = new HashSet<>();
			for (int r=0; r < getNumberOfReactions(); r++)
				if (getProductionStoichiometry(s, r) != 0 || getConsumptionStoichiometry(s, r) != 0)
					ir.add(r);
			involvingReactionsList.add(ImmutableSet.copyOf(ir));
		}
		return involvingReactionsList;
	}

	@Override
	public Set<Integer> getModifiedSpecies(int r) {
		return modifiedSpeciesListSupplier.get().get(r);
	}

	private List<Set<Integer>> computeModifiedSpeciesList() {
		List<Set<Integer>> modifiedSpeciesList = new ArrayList<>(getNumberOfSpecies());
		for (int r=0; r < getNumberOfReactions(); r++) {
			Set<Integer> ms = new HashSet<>();
			for (int s=0; s < getNumberOfSpecies(); s++)
				if (getStoichiometry(s, r) != 0)
					ms.add(s);
			modifiedSpeciesList.add(ImmutableSet.copyOf(ms));
		}
		return modifiedSpeciesList;
	}


	@Override
	public Set<Integer> getModifyingReactions(int s) {
		return modifyingReactionsListSupplier.get().get(s);
	}

	private List<Set<Integer>> computeModifyingReactionsList() {
		List<Set<Integer>> modifyingReactionsList = new ArrayList<>(getNumberOfReactions());
		for (int s=0; s < getNumberOfSpecies(); s++) {
			Set<Integer> mr = new HashSet<>();
			for (int r=0; r < getNumberOfReactions(); r++)
				if (getStoichiometry(s, r) != 0)
					mr.add(r);
			modifyingReactionsList.add(ImmutableSet.copyOf(mr));
		}
		return modifyingReactionsList;
	}
	@Override
	public ReactionNetworkGraph getGraph() {
		return graphSupplier.get();
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
