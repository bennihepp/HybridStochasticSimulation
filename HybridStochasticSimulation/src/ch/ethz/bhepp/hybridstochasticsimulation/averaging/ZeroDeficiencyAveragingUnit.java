package ch.ethz.bhepp.hybridstochasticsimulation.averaging;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.apache.commons.math3.optim.linear.LinearConstraint;
import org.apache.commons.math3.optim.linear.LinearConstraintSet;
import org.apache.commons.math3.optim.linear.LinearObjectiveFunction;
import org.apache.commons.math3.optim.linear.NoFeasibleSolutionException;
import org.apache.commons.math3.optim.linear.Relationship;
import org.apache.commons.math3.optim.linear.SimplexSolver;
import org.apache.commons.math3.optim.nonlinear.scalar.GoalType;
import org.apache.commons.math3.random.RandomDataGenerator;
import org.apache.commons.math3.util.FastMath;
import org.ejml.UtilEjml;
import org.ejml.data.DenseMatrix64F;
import org.ejml.factory.DecompositionFactory;
import org.ejml.factory.SingularValueDecomposition;
import org.ejml.ops.CommonOps;
import org.ejml.ops.SingularOps;
import org.jgrapht.alg.StrongConnectivityInspector;
import org.jgrapht.graph.DirectedSubgraph;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import ch.ethz.bhepp.hybridstochasticsimulation.graphs.ComplexEdge;
import ch.ethz.bhepp.hybridstochasticsimulation.graphs.ComplexGraph;
import ch.ethz.bhepp.hybridstochasticsimulation.graphs.ComplexVertex;
import ch.ethz.bhepp.hybridstochasticsimulation.graphs.SpeciesVertex;
import ch.ethz.bhepp.hybridstochasticsimulation.math.MathUtilities;
import ch.ethz.bhepp.hybridstochasticsimulation.math.distributions.ConstrainedProductFormDistribution;
import ch.ethz.bhepp.hybridstochasticsimulation.math.distributions.IndependentJointDistribution;
import ch.ethz.bhepp.hybridstochasticsimulation.math.distributions.MultivariateDistribution;
import ch.ethz.bhepp.hybridstochasticsimulation.math.distributions.PoissonDistribution;
import ch.ethz.bhepp.hybridstochasticsimulation.math.distributions.UnivariateDistribution;
import ch.ethz.bhepp.hybridstochasticsimulation.models.AdaptiveMSHRNModel;
import ch.ethz.bhepp.hybridstochasticsimulation.models.UnaryBinaryStochasticModel;
import ch.ethz.bhepp.hybridstochasticsimulation.networks.MassActionReactionNetwork;
import ch.ethz.bhepp.hybridstochasticsimulation.networks.ReactionNetworkUtils;

import com.google.common.base.Predicate;

public class ZeroDeficiencyAveragingUnit extends AbstractAveragingUnit implements Serializable {

	private static final long serialVersionUID = 1L;

	private static final Logger logger = LoggerFactory.getLogger(ZeroDeficiencyAveragingUnit.class);

	private RandomDataGenerator rdg;
	private transient List<SubnetworkDescription> zeroDeficiencySubnetworks = null;
	private UnaryBinaryStochasticModel model;
	private boolean sampleFromStationaryDistribution = true;
	private boolean printMessages = false;

	private double[] _currentStationarySolution;
	private Map<SpeciesConservationRelation, MultivariateDistribution> _distributionMap;
	private MultivariateDistribution _unconservedSpeciesDistribution;

	public static ZeroDeficiencyAveragingUnit createCopy(ZeroDeficiencyAveragingUnit averagingUnit, RandomDataGenerator rdg) {
		ZeroDeficiencyAveragingUnit copy = new ZeroDeficiencyAveragingUnit();
		copy.copyFrom(averagingUnit);
		copy.rdg = rdg;
		copy.zeroDeficiencySubnetworks = averagingUnit.zeroDeficiencySubnetworks;
		copy.model = averagingUnit.model;
		copy.printMessages = averagingUnit.printMessages;
		return copy;
	}

	public ZeroDeficiencyAveragingUnit(MassActionReactionNetwork network, Set<Integer> importantSpecies, RandomDataGenerator rdg) {
		super(network, importantSpecies);
		this.model = new UnaryBinaryStochasticModel(network);
		this.rdg = rdg;
	}

	public void setSampleFromStationaryDistribution(boolean sampleFromStationaryDistribution) {
		this.sampleFromStationaryDistribution = sampleFromStationaryDistribution;
	}

	public void setPrintMessages(boolean printMessages) {
		this.printMessages = printMessages;
	}

	protected ZeroDeficiencyAveragingUnit() {
		super();
	}

	@Override
	public void reset() {
		zeroDeficiencySubnetworks = null;
		this.zeroDeficiencySubnetworks = new ArrayList<>();
//		this.zeroDeficiencySubnetworks = findZeroDeficiencySubnetworks();
		super.reset();
	}

	public boolean isZeroDeficiency(SubnetworkDescription subnetworkDescr) {
		MassActionReactionNetwork subnetwork = createSubReactionNetwork(
				network, subnetworkDescr.getSubnetworkSpecies(), subnetworkDescr.getSubnetworkReactions());

		ComplexGraph subnetworkComplexGraph = createComplexGraphOfSubnetwork(subnetwork);
		if (printMessages) {
			StringBuilder sb = new StringBuilder();
			sb.append("Computing deficiency for {");
			for (int s : subnetworkDescr.getSubnetworkSpecies()) {
				SpeciesVertex v = network.getGraph().getSpeciesVertex(s);
				sb.append(v);
				sb.append(", ");
			}
			sb.delete(sb.length() - 2, sb.length());
			sb.append("}");
			logger.info(sb.toString());
		}
		int deficiency = computeDeficiency(subnetwork, subnetworkComplexGraph);
		if (deficiency != 0)
			return false;
		boolean weaklyReversible = isWeaklyReversible(subnetworkComplexGraph);
		if (!weaklyReversible)
			return false;
		boolean irreducible = isIrreducible(subnetwork);
		if (!irreducible)
			return false;

//		subnetworkInformationMap.put(subnetworkInfo, subnetworkInfo);

//		zeroDeficiencySubnetworks.add(subnetworkInfo);

		return true;
	}

	private boolean isIrreducible(MassActionReactionNetwork subnetwork) {
		if (subnetwork.getNumberOfSpecies() == 0)
			return true;
		DenseMatrix64F S = new DenseMatrix64F(subnetwork.getNumberOfSpecies(), subnetwork.getNumberOfReactions());
//		RealMatrix S = new Array2DRowRealMatrix(subnetwork.getNumberOfSpecies(), subnetwork.getNumberOfReactions());
		for (int s=0; s < subnetwork.getNumberOfSpecies(); s++)
			for (int r=0; r < subnetwork.getNumberOfReactions(); r++)
//				S.setEntry(s, r, subnetwork.getStoichiometry(s, r));
				S.set(s, r, subnetwork.getStoichiometry(s, r));
		DenseMatrix64F transS = new DenseMatrix64F(S.numCols, S.numRows);
		CommonOps.transpose(S, transS);
		SingularValueDecomposition<DenseMatrix64F> svd = computeSVD(transS);
		DenseMatrix64F gamma = computeNullSpace(svd);
		boolean[] conservedMask = new boolean[subnetwork.getNumberOfSpecies()];
		for (int c=0; c < gamma.numCols; c++)
			for (int i=0; i < gamma.numRows; i++)
				if (gamma.get(i, c) != 0)
					conservedMask[i] = true;
		int dc = MathUtilities.count(conservedMask);
		if (dc > 0) {
			if (dc == subnetwork.getNumberOfSpecies())
				return true;
			else
				return false;
		} else if (computeRank(S) != subnetwork.getNumberOfSpecies()){
			if (!checkColSpanIntegers(S))
				return false;
			if (!checkColSpanPositiveReals(S))
				return false;
			int[][] nu = subnetwork.getConsumptionStoichiometries();
			int[][] nuPrime = subnetwork.getProductionStoichiometries();
			if (!checkExhaustiveStructure(nu, nuPrime))
				return false;
			if (!checkExhaustiveStructure(nuPrime, nu))
				return false;
			return true;
		}
		return false;
	}

	private boolean checkColSpanPositiveReals(DenseMatrix64F M) {
		// We need to check the feasibility of the following linear program:
		//  { v \in RR^k : M*v = 0 and v_i >= 1 }

		// Define linear objective function.
		double[] qVector = new double[M.numCols];
		for (int i=0; i < M.numCols; i++)
			qVector[i] = 1.0;
		LinearObjectiveFunction objectiveFunction = new LinearObjectiveFunction(qVector, 0.0);

		// Define linear constraints
		LinkedList<LinearConstraint> constraintList = new LinkedList<LinearConstraint>();
		// Constraints: M * v = 0
		for (int i=0; i < M.numRows; i++) {
			for (int k=0; k < M.numCols; k++) {
				qVector[k] = M.get(i, k);
			}
			constraintList.add(new LinearConstraint(qVector, Relationship.EQ, 0.0));
		}
		Arrays.fill(qVector, 0.0);
		// Constraints: v_i >= 1
		for (int k=0; k < M.numRows; k++) {
			qVector[k] = 1.0;
			constraintList.add(new LinearConstraint(qVector, Relationship.GEQ, 1.0));
			qVector[k] = 0.0;
		}
		LinearConstraintSet constraints = new LinearConstraintSet(constraintList);

		// This is probably handled as well from the library
//		InitialGuess ig = new InitialGuess(qVector);

		// Run optimization
		SimplexSolver solver = new SimplexSolver();
		try {
//			PointValuePair pv = solver.optimize(objectiveFunction, constraints, GoalType.MINIMIZE, ig);
			solver.optimize(objectiveFunction, constraints, GoalType.MINIMIZE);
		} catch (NoFeasibleSolutionException | org.apache.commons.math3.exception.TooManyIterationsException e) {
			return false;
		}
		return true;
	}

	private boolean checkColSpanIntegers(DenseMatrix64F M) {
		DenseMatrix64F H = computeHermiteNormalForm(M);
		return H.get(0, 0) != 0;
	}

	private DenseMatrix64F computeHermiteNormalForm(DenseMatrix64F M) {
		int m = M.numRows;
		int n = M.numCols;
		DenseMatrix64F H = new DenseMatrix64F(M);
		// Identity matrix
		DenseMatrix64F U = new DenseMatrix64F(n, n);
		for (int i=0; i < n; i++)
			U.set(i, i, 1);
		int r = FastMath.min(m, n);
		// Location of pivot
		int p = 0;
		for (int i=0; i < r; i++) {
			int nonZeroIndex = -1;
			for (int j=p; j < n; j++)
				if (H.get(i, j) != 0) {
					nonZeroIndex = j;
					break;
				}
			if (nonZeroIndex >= 0) {
				boolean finished = false;
				while (!finished) {
					double minValue = Integer.MAX_VALUE;
					int minValueIndex = -1;
					// Find smallest non-zero entry
					for (int j=p; j < n; j++)
						if (H.get(i, j) != 0 && FastMath.abs(H.get(i, j)) < minValue) {
							minValue = FastMath.abs(H.get(i, j));
							minValueIndex = j;
						}
					if (minValueIndex != p) {
						// exchange columns i and minValueIndex
						swapColumns(H, i, minValueIndex);
						swapColumns(U, i, minValueIndex);
					}
					for (int j=p+1; j < n; j++) {
						int q = (int)FastMath.round(H.get(i, j) / H.get(i, p));
						for (int l=0; l < m; l++) {
							H.set(l, j, H.get(l, j) - q * H.get(l, p));
							U.set(l, j, U.get(l, j) - q * U.get(l, p));
						}
					}
					nonZeroIndex = -1;
					for (int j=p+1; j < n; j++)
						if (H.get(i, j) != 0) {
							nonZeroIndex = j;
							break;
						}
					if (nonZeroIndex < 0) {
						finished = true;
						// Flip sign of ith row if necessary
						if (H.get(i, p) < 0) {
							for (int l=0; l < m; l++) {
								H.set(l, p, - H.get(l, p));
								U.set(l, p, - U.get(l, p));
							}
							// reduce entries to the left of H_{ii}
							for (int j=0; j < p; j++) {
								double q = FastMath.floor(H.get(i, j) / H.get(i, p));
								for (int l=0; l < m; l++) {
									H.set(l, j, H.get(l, j) - q * H.get(l, p));
									U.set(l, j, U.get(l, j) - q * U.get(l, p));
								}
							}
						}
					}
				}
				p++;
			}
		}
		return H;
	}

	private void swapColumns(DenseMatrix64F M, int j1, int j2) {
		for (int i=0; i < M.numRows; i++) {
			double temp = M.get(i, j1);
			M.set(i, j1, M.get(i, j2));
			M.set(i, j2, temp);
		}
	}

	private boolean checkExhaustiveStructure(int[][] nu, int[][] nuPrime) {
		int n = nu.length;
		int m = nu[0].length;
		boolean[] mask = new boolean[n];
		boolean change = true;
		// FIXME: No need to use newMask and it will even be faster.
		while (change) {
			change = false;
			boolean[] newMask = mask.clone();
			for (int k=1; k < m; k++) {
				boolean inMask = true;
				for (int i=0; i < n; i++)
					if (nu[i][k] != 0 && !mask[i]) {
						inMask = false;
					}
				if (inMask) {
					for (int i=0; i < n; i++)
						if (nu[i][k] != 0)
							mask[i] = true;
					change = true;
				}
			}
			mask = newMask;
			boolean allMask = true;
			for (int i=0; i < n; i++)
				if (!mask[i]) {
					allMask = false;
					break;
				}
			if (allMask)
				return true;
		}
		return false;
	}

	public SingularValueDecomposition<DenseMatrix64F> computeSVD(DenseMatrix64F matrix) {
		if (matrix.numRows == 0 || matrix.numCols == 0)
			return null;
		SingularValueDecomposition<DenseMatrix64F> svd = DecompositionFactory.svd(matrix.numRows, matrix.numCols, false, true, false);
		if (!svd.decompose(matrix))
			throw new AveragingException("Unable to perform singular value decomposition");
		return svd;
	}

	public int computeRank(SingularValueDecomposition<DenseMatrix64F> svd) {
		return SingularOps.rank(svd, UtilEjml.EPS);
	}

//	private int computeRank(QRDecomposition qr, double tolerance) {
//		RealMatrix R = qr.getR();
//		int rank = 0;
//		for (int i=0; i < R.getRowDimension(); i++)
//			if (FastMath.abs(R.getEntry(i, i)) > tolerance)
//				rank++;
//		return rank;
//	}

	private DenseMatrix64F computeNullSpace(SingularValueDecomposition<DenseMatrix64F> svd) {
		int rank = SingularOps.rank(svd, UtilEjml.EPS);
		int r = FastMath.min(svd.numRows(), svd.numCols());
		DenseMatrix64F nullSpace = new DenseMatrix64F(svd.numCols(), r - rank);
		return SingularOps.nullSpace(svd, nullSpace, UtilEjml.EPS);
	}

//	private RealMatrix computeNullSpace(QRDecomposition qr, int rank) {
//		RealMatrix Q = qr.getQ();
//		int rows = ;
//		int cols = Q
//		RealMatrix nullSpace = new Array2DRowRealMatrix(M.getRowDimension(), M.getColumnDimension() - rank);
//		for (int c=0; c < nullSpace.getColumnDimension(); c++)
//			nullSpace.setColumnVector(c, Q.getColumnVector(Q.getColumnDimension() - c));
//		return nullSpace;
//	}

	private boolean isWeaklyReversible(ComplexGraph complexGraph) {
		// A chemical reaction network is weakly reversible if and only if all connected components of the complex graph are strongly connected
		boolean weaklyReversible = true;
		List<Set<ComplexVertex>> connectedSets = complexGraph.getConnectedSets();
		for (Set<ComplexVertex> connectedSet : connectedSets) {
			Set<ComplexEdge> connectedSetEdges = new HashSet<>();
			for (ComplexEdge e : complexGraph.edgeSet()) {
				if (connectedSet.contains(complexGraph.getEdgeSource(e)) && connectedSet.contains(complexGraph.getEdgeTarget(e)))
					connectedSetEdges.add(e);
			}
			DirectedSubgraph<ComplexVertex, ComplexEdge> subgraph = new DirectedSubgraph<>(complexGraph, connectedSet, connectedSetEdges);
			StrongConnectivityInspector<ComplexVertex, ComplexEdge> sci = new StrongConnectivityInspector<>(subgraph);
			if (!sci.isStronglyConnected()) {
				weaklyReversible = false;
				break;
			}
		}
		if (printMessages)
			logger.info("  Weakly reversible: " + weaklyReversible);
		return weaklyReversible;
	}

	private int computeDeficiency(MassActionReactionNetwork network, ComplexGraph complexGraph) {
		DenseMatrix64F matrix = ReactionNetworkUtils.createStochiometryMatrix(network);
		int rank = computeRank(matrix);

		List<Set<ComplexVertex>> connectedSets = complexGraph.getConnectedSets();

		int numOfComplexes = complexGraph.vertexSet().size();
		int dimOfStochiometricSubspace = rank;
		int numOfLinkageClasses = connectedSets.size();
		int deficiency = numOfComplexes - numOfLinkageClasses - dimOfStochiometricSubspace;

		if (printMessages) {
			for (Set<ComplexVertex> connectedSet : connectedSets) {
				logger.info("  Connected set:");
				for (ComplexVertex v : connectedSet) {
					logger.info("    " + v);
				}
			}
			logger.info("  Number of complexes: " + numOfComplexes);
			logger.info("  Dimension of stochiometric subspace: " + dimOfStochiometricSubspace);
			logger.info("  Number of linkage classes: " + numOfLinkageClasses);
			logger.info("  Deficiency: " + deficiency);
		}

		return deficiency;
	}

	public int computeRank(DenseMatrix64F matrix) {
		if (matrix.numRows == 0 || matrix.numCols == 0)
			return 0;
		SingularValueDecomposition<DenseMatrix64F> svd = DecompositionFactory.svd(matrix.numRows, matrix.numCols, false, false, false);
		if (!svd.decompose(matrix))
			throw new AveragingException("Unable to perform singular value decomposition");
		return SingularOps.rank(svd, UtilEjml.EPS);
	}

	private ComplexGraph createComplexGraphOfSubnetwork(MassActionReactionNetwork subnetwork) {
		ComplexGraph complexGraph = ComplexGraph.createFromReactionNetwork(subnetwork);
		return complexGraph;
	}

	@Override
	public void sampleSubnetworkState(double t, double[] x, SubnetworkInformation subnetworkInfo) {
		double[] stationarySolution;
		if (_currentStationarySolution != null)
			stationarySolution = _currentStationarySolution;
		else
			stationarySolution = subnetworkInfo.computeStationarySolution(t, x);

		if (sampleFromStationaryDistribution) {

			for (SpeciesConservationRelation relation : subnetworkInfo.getConservedSpeciesRelations()) {
				// Sample from the stationary distribution for each independent reducible subnetwork
				MultivariateDistribution conservedSpeciesDistribution;
				if (_distributionMap != null)
					conservedSpeciesDistribution = _distributionMap.get(relation);
				else
					conservedSpeciesDistribution = computeConservedSpeciesDistribution(subnetworkInfo, relation, stationarySolution);
				double[] y = conservedSpeciesDistribution.sample();
				List<Integer> speciesList = relation.getConservedSpeciesList();
				for (int i=0; i < speciesList.size(); i++) {
					int s = speciesList.get(i);
					x[s] = y[i];
				}
			}
	
			MultivariateDistribution unconservedSpeciesDistribution;
			if (_unconservedSpeciesDistribution != null)
				unconservedSpeciesDistribution = _unconservedSpeciesDistribution;
			else
				unconservedSpeciesDistribution = computeUnconservedSpeciesDistribution(subnetworkInfo, stationarySolution);
			double[] y = unconservedSpeciesDistribution.sample();
			List<Integer> unconservedSpecies = subnetworkInfo.getUnconservedSpecies();
			for (int i=0; i < subnetworkInfo.getUnconservedSpecies().size(); i++) {
				int s = unconservedSpecies.get(i);
				x[s] = y[i];
			}

		} else {

			int i = 0;
//			for (int i=0; i < subnetworkInfo.getSubnetworkSpecies().size(); i++) {
			Iterator<Integer> it = subnetworkInfo.getSubnetworkSpecies().iterator();
			while (it.hasNext()) {
				int j = it.next();
//				x[subnetworkInfo.getIndexMap().get(i)] = stationarySolution[i];
				x[j] = stationarySolution[i];
				++i;
			}

		}
	}

	@Override
	protected double[] computeStationaryFirstMoments(double t, double[] x, SubnetworkInformation subnetworkInfo) {
		double[] firstMoments = new double[subnetworkInfo.getSubnetworkSpecies().size()];

		for (SpeciesConservationRelation relation : subnetworkInfo.getConservedSpeciesRelations()) {
			// Sample from the stationary distribution for each independent reducible subnetwork
			MultivariateDistribution distribution = _distributionMap.get(relation);
			double[] fm = distribution.getFirstMoment();
			List<Integer> speciesList = relation.getConservedSpeciesList();
			for (int i=0; i < speciesList.size(); i++) {
//				int s = speciesList.get(i);
				if (Double.isNaN(fm[i]))
					throw new RuntimeException();
				firstMoments[i] = fm[i];
			}
		}

		double[] fm = _unconservedSpeciesDistribution.getFirstMoment();
		List<Integer> unconservedSpecies = subnetworkInfo.getUnconservedSpecies();
		for (int i=0; i < subnetworkInfo.getUnconservedSpecies().size(); i++) {
			int s = unconservedSpecies.get(i);
			firstMoments[s] = fm[i];
		}

		return firstMoments;
	}

	@Override
	protected double[][] computeStationarySecondMoments(double t, double[] x, SubnetworkInformation subnetworkInfo) {
		double[][] secondMoments = new double[subnetworkInfo.getSubnetworkSpecies().size()][subnetworkInfo.getSubnetworkSpecies().size()];

		for (SpeciesConservationRelation relation : subnetworkInfo.getConservedSpeciesRelations()) {
			// Sample from the stationary distribution for each independent reducible subnetwork
			MultivariateDistribution distribution = _distributionMap.get(relation);
			double[][] sm = distribution.getSecondMoment();
			List<Integer> speciesList = relation.getConservedSpeciesList();
			for (int i=0; i < speciesList.size(); i++) {
				for (int j=0; j < speciesList.size(); j++) {
//					int s1 = speciesList.get(i);
//					int s2 = speciesList.get(j);
					secondMoments[i][j] = sm[i][j];
				}
			}
		}

		double[][] sm = _unconservedSpeciesDistribution.getSecondMoment();
		List<Integer> unconservedSpecies = subnetworkInfo.getUnconservedSpecies();
		for (int i=0; i < subnetworkInfo.getUnconservedSpecies().size(); i++) {
			for (int j=0; j < subnetworkInfo.getUnconservedSpecies().size(); j++) {
				int s1 = unconservedSpecies.get(i);
				int s2 = unconservedSpecies.get(j);
				secondMoments[s1][s2] = sm[i][j];
			}
		}

		// Compute all other second moments

		return secondMoments;
	}

	@Override
	public Predicate<SubnetworkDescription> getSubnetworkFilter() {
		return new Predicate<SubnetworkDescription>() {

			@Override
			public boolean apply(SubnetworkDescription subnetworkDescr) {
				return isZeroDeficiency(subnetworkDescr);
			}

		};
	}

	@Override
	public void updateSubnetworkState(AdaptiveMSHRNModel model, double t, double[] x, SubnetworkDescription subnetwork) {
		SubnetworkInformation subnetworkInfo = subnetworkInformationMap.get(subnetwork);

		_currentStationarySolution = subnetworkInfo.computeStationarySolution(t, x);
		_distributionMap = new HashMap<>(subnetworkInfo.getConservedSpeciesRelations().size());
		for (SpeciesConservationRelation relation : subnetworkInfo.getConservedSpeciesRelations()) {
			MultivariateDistribution distribution = computeConservedSpeciesDistribution(subnetworkInfo, relation, _currentStationarySolution);
			_distributionMap.put(relation, distribution);
		}
		_unconservedSpeciesDistribution = computeUnconservedSpeciesDistribution(subnetworkInfo, _currentStationarySolution);

		super.updateSubnetworkState(model, t, x, subnetworkInfo);
		_unconservedSpeciesDistribution = null;
		_distributionMap = null;
		_currentStationarySolution = null;
	}

	private MultivariateDistribution computeUnconservedSpeciesDistribution(SubnetworkInformation subnetworkInfo, double[] stationarySolution) {
		List<UnivariateDistribution> poissonDistributions = new ArrayList<>(subnetworkInfo.getUnconservedSpecies().size());
		for (int i : subnetworkInfo.getUnconservedSpecies()) {
			double lambda = stationarySolution[i];
			poissonDistributions.add(new PoissonDistribution(rdg, lambda));
		}
		MultivariateDistribution distribution = new IndependentJointDistribution(poissonDistributions);
		return distribution;
	}

	private MultivariateDistribution computeConservedSpeciesDistribution(
			SubnetworkInformation subnetworkInfo,
			SpeciesConservationRelation relation, double[] stationarySolution) {
		// Sample from the stationary distribution for each independent reducible subnetwork
		Map<Integer, Integer> indexMap = subnetworkInfo.getIndexMap();
		List<Integer> speciesList = relation.getConservedSpeciesList();
		DenseMatrix64F lcVector = relation.getLinearCombination();

		double nDouble = 0.0;
		int[] alpha = new int[speciesList.size()];
		double[] c = new double[speciesList.size()];
		for (int j=0; j < speciesList.size(); j++) {
			int s = speciesList.get(j);
			int k = indexMap.get(s);
			alpha[j] = (int)FastMath.round(lcVector.get(j, 0));
			c[j] = stationarySolution[k];
			nDouble += alpha[j] * c[j];
		}

		int n = (int)FastMath.round(nDouble);

//		if (n == 0)
//			// All conserved species are zero so we don't need to sample anything
//			continue;

		return new ConstrainedProductFormDistribution(rdg, n, c, alpha);
	}

}
