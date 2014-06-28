/**
 * @author      Benjamin Hepp <benjamin.hepp@bsse.ethz.ch>
 */
package ch.ethz.bhepp.hybridstochasticsimulation.averaging;

import static com.google.common.base.Preconditions.checkArgument;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.apache.commons.math3.util.FastMath;
import org.ejml.UtilEjml;
import org.ejml.data.DenseMatrix64F;
import org.ejml.factory.DecompositionFactory;
import org.ejml.factory.SingularValueDecomposition;
import org.ejml.ops.CommonOps;
import org.ejml.ops.SingularOps;

import cern.colt.matrix.DoubleMatrix1D;
import cern.colt.matrix.DoubleMatrix2D;
import cern.colt.matrix.impl.DenseDoubleMatrix1D;
import cern.colt.matrix.impl.SparseDoubleMatrix2D;
import cern.colt.matrix.linalg.Algebra;
import cern.jet.math.Functions;
import ch.ethz.bhepp.hybridstochasticsimulation.graphs.SpeciesVertex;
import ch.ethz.bhepp.hybridstochasticsimulation.math.MathUtilities;
import ch.ethz.bhepp.hybridstochasticsimulation.models.UnaryBinaryStochasticModel;
import ch.ethz.bhepp.hybridstochasticsimulation.networks.MassActionReactionNetwork;

import com.google.common.base.Predicate;

public class FiniteMarkovChainAveragingUnit extends AbstractAveragingUnit {

	// TODO: StateEnumerator could be implemented more efficiently
	// TODO: Reuse StateEnumerator and transitionMatrix etc. instead of computing again and again

	public static FiniteMarkovChainAveragingUnit createCopy(FiniteMarkovChainAveragingUnit averagingUnit) {
		FiniteMarkovChainAveragingUnit copy = new FiniteMarkovChainAveragingUnit();
		copy.copyFrom(averagingUnit);
		copy.finiteMarkovChainSubnetworks = averagingUnit.finiteMarkovChainSubnetworks;
		copy.subnetworkInformationMap = averagingUnit.subnetworkInformationMap;
		copy._stopIfAveragingBecomesInvalid = averagingUnit._stopIfAveragingBecomesInvalid;
		copy._warnIfAveragingBecomesInvalid = averagingUnit._warnIfAveragingBecomesInvalid;
		return copy;
	}

	private List<Set<SpeciesVertex>> finiteMarkovChainSubnetworks = null;
	private Map<Set<SpeciesVertex>, SubnetworkInformation> subnetworkInformationMap;
	// TODO: Cache the transition matrix?
//	private Map<Set<SpeciesVertex>, Matrix> averagingCandidateTransitionMatrixMap;
	private boolean _stopIfAveragingBecomesInvalid = true;
	private boolean _warnIfAveragingBecomesInvalid = true;
//	private double zeroEigenvalueTolerance = 1e-6;
	private double zeroSumTolerance = 1e-6;

	public FiniteMarkovChainAveragingUnit(MassActionReactionNetwork network, Set<Integer> importantSpecies) {
		super(network, importantSpecies);
		this.subnetworkInformationMap = new HashMap<>();
	}

	protected FiniteMarkovChainAveragingUnit() {
		super();
	}

	public void stopIfAveragingBecomesInvalid(boolean stop) {
		_stopIfAveragingBecomesInvalid = stop;
	}

	public void warnIfAveragingBecomesInvalid(boolean warn) {
		_warnIfAveragingBecomesInvalid = warn;
	}

	public void setZeroSumTolerance(double tol) {
		this.zeroSumTolerance = tol;
	}

//	public void setZeroEigenvalueTolerance(double tol) {
//		this.zeroEigenvalueTolerance = tol;
//	}

//	private List<Set<SpeciesVertex>> findFiniteMarkovChainSubnetworks() {
//		List<Set<SpeciesVertex>> finiteMarkovChainSubnetworks = new ArrayList<Set<SpeciesVertex>>();
//		for (Set<SpeciesVertex> subnetworkSpecies : enumerateSubnetworks()) {
//			MassActionReactionNetwork subnetwork = createSubReactionNetwork(network, subnetworkSpecies);
//
//			List<SpeciesConservationRelation> conservedSpeciesRelations = getConservedSpeciesRelations(subnetwork);
//			// TODO: for now we only allow subnetworks that are fully related by a single conservation relation
//			if (conservedSpeciesRelations.size() != 1)
//				continue;
////			Set<SpeciesVertex> unconservedSpeciesSet = new HashSet<>(subnetworkSpecies);
////			for (SpeciesConservationRelation relation : conservedSpeciesRelations)
////				unconservedSpeciesSet.removeAll(relation.getConservedSpeciesList());
//////			// Check if this a Finite Markov Chain (i.e. all species in the subnetwork are involved in a conservation relation)
////			if (unconservedSpeciesSet.size() > 0)
////				continue;
//			SpeciesConservationRelation conservedSpeciesRelation = conservedSpeciesRelations.get(0);
//
//			UnaryBinaryStochasticModel subnetworkModel = new UnaryBinaryStochasticModel(subnetwork);
//
//			SubnetworkInformation subnetworkInfo = new SubnetworkInformation();
//			subnetworkInfo.setConservedSpeciesRelation(conservedSpeciesRelation);
//			subnetworkInfo.setModel(subnetworkModel);
//			Map<SpeciesVertex, Integer> subnetworkIndexMap = new HashMap<>();
//			for (int i=0; i < subnetwork.getNumberOfSpecies(); i++) {
//				SpeciesVertex v = subnetwork.getGraph().getSpeciesVertex(i);
//				subnetworkIndexMap.put(v, i);
//			}
//			subnetworkInfo.setIndexMap(subnetworkIndexMap);
//			subnetworkInformationMap.put(subnetworkSpecies, subnetworkInfo);
//
//			finiteMarkovChainSubnetworks.add(subnetworkSpecies);
//		}
//		return finiteMarkovChainSubnetworks;
//	}

//	private MassActionReactionNetwork createSubReactionNetwork(MassActionReactionNetwork network, Set<SpeciesVertex> subnetworkSpecies) {
//		Set<Integer> subSpeciesIndices = new HashSet<Integer>(subnetworkSpecies.size());
//		for (SpeciesVertex v : subnetworkSpecies)
//			subSpeciesIndices.add(v.getSpecies());
//		DefaultUnaryBinaryReactionNetwork subnetwork = DefaultUnaryBinaryReactionNetwork.createSubnetwork(network, subSpeciesIndices);
//		return subnetwork;
//	}

//	@Override
//	public List<Set<SpeciesVertex>> findAveragingCandidates(double t, double[] x, Predicate<Set<SpeciesVertex>> filter) {
//		if (finiteMarkovChainSubnetworks == null)
//			this.finiteMarkovChainSubnetworks = findFiniteMarkovChainSubnetworks();
////		averagingCandidateTransitionMatrixMap = new HashMap<>();
//		// First find a list of subnetworks that could be averaged
//		List<Set<SpeciesVertex>> averagingCandidates = new ArrayList<Set<SpeciesVertex>>();
//		for (Set<SpeciesVertex> subnetworkSpecies : finiteMarkovChainSubnetworks) {
//			if (filter.apply(subnetworkSpecies)) {
//				SubnetworkInformation subnetworkInfo = subnetworkInformationMap.get(subnetworkSpecies);
//				// Colt
//				DoubleMatrix2D transitionMatrix = computeTransitionMatrix(t, x, subnetworkInfo.getConservedSpeciesRelation(), subnetworkInfo.getModel());
//				// La4j
////				Matrix transitionMatrix = computeTransitionMatrix(t, x, subnetworkInfo.getConservedSpeciesRelation(), subnetworkInfo.getModel());
//				if (isIrreducible(transitionMatrix)) {
////					averagingCandidateTransitionMatrixMap.put(subnetworkSpecies, transitionMatrix);
//					averagingCandidates.add(subnetworkSpecies);
//				}
//			}
//		}
//		return averagingCandidates;
//	}

//	@Override
//	public void resampleFromStationaryDistribution(double t, double[] x, Set<SpeciesVertex> subnetworkSpecies) {
//		SubnetworkInformation subnetworkInfo = subnetworkInformationMap.get(subnetworkSpecies);
//		SpeciesConservationRelation relation = subnetworkInfo.getConservedSpeciesRelation();
//		UnaryBinaryStochasticModel subnetworkModel = subnetworkInfo.getModel();
//		DoubleMatrix2D transitionMatrix = computeTransitionMatrix(t, x, relation, subnetworkModel);
//		DoubleMatrix1D stationaryDistribution = computeStationaryDistribution(transitionMatrix);
////		Matrix transitionMatrix = computeTransitionMatrix(t, x, relation, subnetworkModel);
////		Vector stationaryDistribution = computeStationaryDistribution(transitionMatrix);
//		int stateIndex = RandomDataUtilities.sampleFromProbabilityMassFunction(rdg, stationaryDistribution.toArray());
////		double[] stationaryDistributionArray = convertVectorToArray(stationaryDistribution);
////		int stateIndex = RandomDataUtilities.sampleFromProbabilityMassFunction(rdg, stationaryDistributionArray);
//		SpeciesConservationRelationStateEnumerator stateEnumerator = new SpeciesConservationRelationStateEnumerator(x, relation);
//		Iterator<double[]> it = stateEnumerator.iterator();
//		for (int i=0; i < stateIndex; i++)
//			it.next();
//		double[] state = it.next();
//		Map<SpeciesVertex, Integer> indexMap = subnetworkInfo.getIndexMap();
//		for (SpeciesVertex v : subnetworkSpecies) {
//			x[v.getSpecies()] = state[indexMap.get(v)];
//		}
//	}

//	private double[] convertVectorToArray(Vector stationaryDistribution) {
//		double[] array = new double[stationaryDistribution.length()];
//		for (int i=0; i < array.length; i++)
//			array[i] = stationaryDistribution.get(i);
//		return array;
//	}

	private List<SpeciesConservationRelation> getConservedSpeciesRelations(MassActionReactionNetwork network) {
		DenseMatrix64F nullSpace = computeNullSpaceOfReactionNetwork(network);
		List<SpeciesConservationRelation> conservedSpeciesRelations = new ArrayList<>();
		for (int col=0; col < nullSpace.numCols; col++) {
			List<Integer> speciesList = new ArrayList<>();
			boolean noConservationRelation = false;
			double sign = 0.0;
			for (int row=0; row < nullSpace.numRows; row++) {
				double v = nullSpace.get(row, col);
				if (sign == 0.0 && v != 0.0)
					sign = v > 0.0 ? 1.0 : -1.0;
				if (v * sign < 0.0) {
					noConservationRelation = true;
					break;
				}
				if (v != 0.0)
					speciesList.add(row);
			}
			if (noConservationRelation)
				continue;
			DenseMatrix64F lcVector = CommonOps.extract(nullSpace, 0, nullSpace.numRows, col, col + 1);
			CommonOps.scale(sign, lcVector);
			double elementMinNeqZero = Double.MAX_VALUE;
			for (int i=0; i < lcVector.numRows; i++) {
				double element = lcVector.get(i);
				if (element < elementMinNeqZero && element > 0.0)
					elementMinNeqZero = element;
			}
			CommonOps.divide(elementMinNeqZero, lcVector);
			conservedSpeciesRelations.add(new SpeciesConservationRelation(speciesList, lcVector));
		}
		return conservedSpeciesRelations;
	}

	private DenseMatrix64F computeNullSpaceOfReactionNetwork(MassActionReactionNetwork network) {
		DenseMatrix64F matrix = createStochiometryMatrix(network);
		if (matrix.numRows == 0 || matrix.numCols == 0)
			return new DenseMatrix64F(0, 0);
		SingularValueDecomposition<DenseMatrix64F> svd = DecompositionFactory.svd(matrix.numRows, matrix.numCols, false, true, false);
		if (!svd.decompose(matrix))
			throw new AveragingException("Unable to perform singular value decomposition");
		DenseMatrix64F nullSpace = new DenseMatrix64F(matrix.numRows, matrix.numCols);
		SingularOps.nullSpace(svd, nullSpace, UtilEjml.EPS);
		return nullSpace;
	}

	private DenseMatrix64F createStochiometryMatrix(MassActionReactionNetwork network) {
		DenseMatrix64F matrix = new DenseMatrix64F(network.getNumberOfReactions(), network.getNumberOfSpecies());
		for (int r=0; r < network.getNumberOfReactions(); r++) {
			int[] productionStochtiometries = network.getProductionStoichiometries(r);
			int[] consumptionStochtiometries = network.getConsumptionStoichiometries(r);
			for (int s=0; s < network.getNumberOfSpecies(); s++) {
				int stochiometry = productionStochtiometries[s] - consumptionStochtiometries[s];
				matrix.set(r, s, stochiometry);
			}
		}
		return matrix;
	}

//	@Override
//	protected void computeAverageStationaryStateOfSubnetworks(double t, double[] x, List<Set<SpeciesVertex>> subnetworksToAverage) {
//		// TODO: Perform correct computation of average state
//		for (Set<SpeciesVertex> subnetworkSpecies : subnetworksToAverage) {
//			SubnetworkInformation subnetworkInfo = subnetworkInformationMap.get(subnetworkSpecies);
//			SpeciesConservationRelation conservedSpeciesRelation = subnetworkInfo.getConservedSpeciesRelation();
//			UnaryBinaryStochasticModel subnetworkModel = subnetworkInfo.getModel();
//			Map<SpeciesVertex, Integer> indexMap = subnetworkInfo.getIndexMap();
//			double[] subXAverage = computeAverageStateOfSubnetwork(t, x, conservedSpeciesRelation, subnetworkModel);
//			for (SpeciesVertex v : conservedSpeciesRelation.getConservedSpeciesList()) {
//				x[v.getSpecies()] = subXAverage[indexMap.get(v)];
//			}
//		}
//	}

//	private double[] computeAverageStateOfSubnetwork(double t, double[] x, SpeciesConservationRelation relation,
//			UnaryBinaryStochasticModel subnetworkModel) {
//		DoubleMatrix2D transitionMatrix = computeTransitionMatrix(t, x, relation, subnetworkModel);
////		Matrix transitionMatrix = computeTransitionMatrix(t, x, relation, subnetworkModel);
//		EigenvalueDecomposition eigenDecomposition = new EigenvalueDecomposition(transitionMatrix);
////		MatrixDecompositor eigenDecompositor = transitionMatrix.withDecompositor(DecompositorFactory.EIGEN);
////		Matrix[] pdp = eigenDecompositor.decompose(LinearAlgebra.CRS_FACTORY);
////		Matrix P = pdp[0];
////		Matrix D = pdp[1];  
//		// TODO: Is this even possible?
//		boolean irreducible = isIrreducible(eigenDecomposition);
////		boolean irreducible = isIrreducible(P, D);
//		if (!irreducible) {
//			String msg = "Averaging of finite markov chain subnetwork became invalid. The finite markov chain is not irreducible anymore. at t=" + t;
//			if (_stopIfAveragingBecomesInvalid)
//				throw new AveragingException(msg);
//			if (_warnIfAveragingBecomesInvalid)
//				logger.warn(msg);
//		}
//		DoubleMatrix1D stationaryDistribution = computeStationaryDistribution(eigenDecomposition);
////		Vector stationaryDistribution = computeStationaryDistribution(P, D);
//		double[] xAverage = computeStationaryAverage(x, relation, stationaryDistribution);
//		return xAverage;
//	}

//	private boolean isIrreducible(DoubleMatrix2D transitionMatrix) {
////	private boolean isIrreducible(Matrix transitionMatrix) {
//		EigenvalueDecomposition eigenDecomposition = new EigenvalueDecomposition(transitionMatrix);
////		double det = transitionMatrix.determinant();
////		LinearSystemSolver solver = transitionMatrix.withSolver(LinearAlgebra.LEAST_SQUARES);
////		Vector b = new BasicVector(transitionMatrix.rows());
////		Vector q = solver.solve(b);
////		MatrixDecompositor eigenDecompositor = transitionMatrix.withDecompositor(DecompositorFactory.EIGEN);
////		Matrix[] pdp = eigenDecompositor.decompose(LinearAlgebra.CRS_FACTORY);
////		Matrix P = pdp[0];
////		Matrix D = pdp[1]; 
//		return isIrreducible(eigenDecomposition);
////		return isIrreducible(P, D);
//	}

//	private boolean isIrreducible(EigenvalueDecomposition eigenDecomposition) {
////	private boolean isIrreducible(Matrix P, Matrix D) {
//		// To check for irreducibility we check whether the dimension of the nullspace is <= 1
//		DoubleMatrix2D D = eigenDecomposition.getD();
//		int dimensionOfNullspace = 0;
//		for (int i=0; i < D.rows(); i++)
//			if (FastMath.abs(D.get(i, i)) <= zeroEigenvalueTolerance)
//				dimensionOfNullspace++;
//		boolean irreducible = dimensionOfNullspace <= 1;
//		return irreducible;
//	}

	private DoubleMatrix2D computeTransitionMatrix(double t, double[] x, SpeciesConservationRelation relation, UnaryBinaryStochasticModel subnetworkModel) {
//	private Matrix computeTransitionMatrix(double t, double[] x, SpeciesConservationRelation relation, UnaryBinaryStochasticModel subnetworkModel) {
		// TODO: check for species numbers > 2
		SpeciesConservationRelationStateEnumerator stateEnumerator = new SpeciesConservationRelationStateEnumerator(x, relation);
		int numOfStates = stateEnumerator.getNumberOfStates();
		DoubleMatrix2D transitionMatrix = new SparseDoubleMatrix2D(numOfStates, numOfStates);
//		double[] diagonal = new double[numOfStates];
//		double[] aboveDiagonal = new double[numOfStates - 1];
//		double[] belowDiagonal = new double[numOfStates - 1];
//		Matrix transitionMatrix = new CRSMatrix(numOfStates, numOfStates);
//		Matrix transitionMatrix = crs.createMatrix(numOfStates, numOfStates);
		int stateIndex = 0;
		double[] propensities = new double[subnetworkModel.getNumberOfReactions()];
		for (double[] state : stateEnumerator) {
			subnetworkModel.computePropensities(stateIndex, state, propensities);
			double propensity = MathUtilities.sum(propensities);
			transitionMatrix.set(stateIndex, stateIndex, -propensity);
			for (int r=0; r < subnetworkModel.getNumberOfReactions(); r++) {
				double newPropensity = propensities[r];
				if (newPropensity > 0) {
					double[] newState = state.clone();
					subnetworkModel.changeState(r, t, newState);
					int newStateIndex = stateEnumerator.getStateIndex(newState);
//					transitionMatrix.setQuick(newStateIndex, stateIndex, newPropensity);
					transitionMatrix.set(newStateIndex, stateIndex, newPropensity);
				}
			}
			stateIndex++;
		}
		return transitionMatrix;
	}

	private DoubleMatrix1D computeStationaryDistribution(DoubleMatrix2D transitionMatrix) {
		Algebra algebra = new Algebra();
		DoubleMatrix2D A = algebra.subMatrix(transitionMatrix, 0, transitionMatrix.rows() - 1, 0, transitionMatrix.columns() - 2);
//		DoubleMatrix2D fullRankMatrix = transitionMatrix.viewPart(0, 0, transitionMatrix.rows(), transitionMatrix.columns() - 1);
		DoubleMatrix2D b = algebra.subMatrix(transitionMatrix, 0, transitionMatrix.rows() - 1, transitionMatrix.columns() - 1, transitionMatrix.columns() - 1);
		DoubleMatrix2D x = algebra.solve(A, b);
		DoubleMatrix1D solution = new DenseDoubleMatrix1D(transitionMatrix.rows());
		for (int i=0; i < x.rows(); i++)
			solution.set(i, x.get(i, 0));
		solution.set(solution.size() - 1, 1.0);
		double sum = solution.aggregate(Functions.plus, Functions.abs);
		solution = solution.assign(Functions.div(sum));
		// Consistency check
		sum = solution.aggregate(Functions.plus, Functions.abs);
		checkArgument(FastMath.abs(sum - 1.0) <= zeroSumTolerance);
		sum = solution.aggregate(Functions.plus, Functions.identity);
		if (sum < 0)
			solution = solution.assign(Functions.mult(-1));
		return solution;
	}

//		private DoubleMatrix1D computeStationaryDistribution(DoubleMatrix2D transitionMatrix) {
////		private Vector computeStationaryDistribution(Matrix transitionMatrix) {
//			EigenvalueDecomposition eigenDecomposition = new EigenvalueDecomposition(transitionMatrix);
////			MatrixDecompositor eigenDecompositor = transitionMatrix.withDecompositor(DecompositorFactory.EIGEN);
////			Matrix[] pdp = eigenDecompositor.decompose(LinearAlgebra.CRS_FACTORY);
////			Matrix P = pdp[0];
////			Matrix D = pdp[1];
//			return computeStationaryDistribution(eigenDecomposition);
////			return computeStationaryDistribution(P, D);
//		}

//	private DoubleMatrix1D computeStationaryDistribution(EigenvalueDecomposition eigenDecomposition) {
////	private Vector computeStationaryDistribution(Matrix P, Matrix D) {
//		DoubleMatrix2D D = eigenDecomposition.getD();
//		DoubleMatrix2D P = eigenDecomposition.getV();
//		// Find eigenvector with eigenvalue 0
//		int zeroEigenvalueIndex = -1;
//		for (int i=0; i < D.rows(); i++)
//			if (FastMath.abs(D.get(i, i)) <= zeroEigenvalueTolerance) {
//				zeroEigenvalueIndex = i;
//				break;
//			}
//		DoubleMatrix1D zeroEigenVector = P.viewColumn(zeroEigenvalueIndex);
//		double sum = zeroEigenVector.aggregate(Functions.plus, Functions.identity);
//		DoubleMatrix1D result = zeroEigenVector.copy();
//		result = result.assign(Functions.div(sum));
////		Vector zeroEigenVector = P.getColumn(zeroEigenvalueIndex);
////		double sum = zeroEigenVector.sum();
////		Vector result = zeroEigenVector.divide(sum);
//		// Consistency check
//		sum = result.aggregate(Functions.plus, Functions.abs);
////		sum = result.sum();
//		checkArgument(FastMath.abs(sum - 1.0) <= zeroSumTolerance);
//		return result;
//	}

	@Override
	public void sampleSubnetworkState(double t, double[] x, SubnetworkDescription subnetwork) {
		// TODO Auto-generated method stub
		
	}

	@Override
	public Predicate<SubnetworkDescription> getSubnetworkFilter() {
		return new Predicate<SubnetworkDescription>() {

			@Override
			public boolean apply(SubnetworkDescription subnetworkDescr) {
				return hasFiniteStateSpace(subnetworkDescr);
			}

		};
	}

	private boolean hasFiniteStateSpace(SubnetworkDescription subnetworkDescr) {
		MassActionReactionNetwork subnetwork = createSubReactionNetwork(
				network, subnetworkDescr.getSubnetworkSpecies(), subnetworkDescr.getSubnetworkReactions());

		List<SpeciesConservationRelation> conservedSpeciesRelations = getConservedSpeciesRelations(subnetwork);

		// TODO: for now we only allow subnetworks that are fully related by a single conservation relation
		return conservedSpeciesRelations.size() == 1;
	}

	private DoubleMatrix1D computeStationaryDistribution(double t, double[] x, SubnetworkInformation subnetworkInfo) {
		SpeciesConservationRelation conservedSpeciesRelation = subnetworkInfo.getConservedSpeciesRelations().get(0);
		UnaryBinaryStochasticModel stochasticModel = subnetworkInfo.getStochasticModel();
		DoubleMatrix2D transitionMatrix = computeTransitionMatrix(t, x, conservedSpeciesRelation, stochasticModel);
//		Matrix transitionMatrix = computeTransitionMatrix(t, x, relation, subnetworkModel);

//		MLDouble transitionMatrixML = MatlabDataExporter.buildDouble("transitionMatrix", transitionMatrix.toArray());
//		SpeciesConservationRelationStateEnumerator stateEnumerator = new SpeciesConservationRelationStateEnumerator(x, conservedSpeciesRelation);
//		double[][] states = new double[stateEnumerator.getNumberOfStates()][stochasticModel.getNumberOfSpecies()];
//		int i = 0;
//		for (double[] state : stateEnumerator) {
//			states[i] = state.clone();
//			i++;
//		}
//		MLDouble statesML = MatlabDataExporter.buildDouble("states", states);
//		ImmutableList<MLArray> matlabData = ImmutableList.<MLArray>of(transitionMatrixML, statesML);
//		try {
//			MatlabDataExporter.writeMatlabDataToFile(new File("/home/bhepp/test.mat"), matlabData);
//		} catch (IOException e) {
//			// TODO Auto-generated catch block
//			e.printStackTrace();
//		}

		DoubleMatrix1D stationaryDistribution = computeStationaryDistribution(transitionMatrix);

//		EigenvalueDecomposition eigenDecomposition = new EigenvalueDecomposition(transitionMatrix);
////		MatrixDecompositor eigenDecompositor = transitionMatrix.withDecompositor(DecompositorFactory.EIGEN);
////		Matrix[] pdp = eigenDecompositor.decompose(LinearAlgebra.CRS_FACTORY);
////		Matrix P = pdp[0];
////		Matrix D = pdp[1];  
//		// TODO: Is this even possible?
//		boolean irreducible = isIrreducible(eigenDecomposition);
////		boolean irreducible = isIrreducible(P, D);
//		if (!irreducible) {
//			String msg = "Averaging of finite markov chain subnetwork became invalid. The finite markov chain is not irreducible anymore. at t=" + t;
//			if (_stopIfAveragingBecomesInvalid)
//				throw new AveragingException(msg);
//			if (_warnIfAveragingBecomesInvalid)
//				logger.warn(msg);
//		}
//		DoubleMatrix1D stationaryDistribution = computeStationaryDistribution(eigenDecomposition);
//		Vector stationaryDistribution = computeStationaryDistribution(P, D);

		return stationaryDistribution;
	}

	@Override
	protected double[] computeStationaryFirstMoments(double t, double[] x, SubnetworkInformation subnetworkInfo) {
		SpeciesConservationRelation conservedSpeciesRelation = subnetworkInfo.getConservedSpeciesRelations().get(0);
		DoubleMatrix1D stationaryDistribution = computeStationaryDistribution(t, x, subnetworkInfo);
		double[] firstMoments = computeStationaryFirstMoments(x, conservedSpeciesRelation, stationaryDistribution);
		return firstMoments;
	}

	private double[] computeStationaryFirstMoments(double[] x, SpeciesConservationRelation relation, DoubleMatrix1D stationaryDistribution) {
//	private double[] computeStationaryAverage(double[] x, SpeciesConservationRelation relation, Vector stationaryDistribution) {
		SpeciesConservationRelationStateEnumerator stateEnumerator = new SpeciesConservationRelationStateEnumerator(x, relation);
		double[] firstMoments = new double[relation.getConservedSpeciesList().size()];
		int stateIndex = 0;
		for (double[] state : stateEnumerator) {
			double p = stationaryDistribution.get(stateIndex);
			for (int i=0; i < firstMoments.length; i++)
				firstMoments[i] += p * state[i];
			stateIndex++;
		}
		return firstMoments;
	}

	@Override
	protected double[][] computeStationarySecondMoments(double t, double[] x, SubnetworkInformation subnetworkInfo) {
		SpeciesConservationRelation conservedSpeciesRelation = subnetworkInfo.getConservedSpeciesRelations().get(0);
		DoubleMatrix1D stationaryDistribution = computeStationaryDistribution(t, x, subnetworkInfo);
		double[][] secondMoments = computeStationarySecondMoments(x, conservedSpeciesRelation, stationaryDistribution);
		return secondMoments;
	}

	private double[][] computeStationarySecondMoments(double[] x, SpeciesConservationRelation relation, DoubleMatrix1D stationaryDistribution) {
//	private double[] computeStationaryAverage(double[] x, SpeciesConservationRelation relation, Vector stationaryDistribution) {
		SpeciesConservationRelationStateEnumerator stateEnumerator = new SpeciesConservationRelationStateEnumerator(x, relation);
		double[][] secondMoments = new double[relation.getConservedSpeciesList().size()][relation.getConservedSpeciesList().size()];
		int stateIndex = 0;
		for (double[] state : stateEnumerator) {
			double p = stationaryDistribution.get(stateIndex);
			for (int i=0; i < secondMoments.length; i++)
				for (int j=0; j < secondMoments[j].length; j++)
					secondMoments[i][j] += p * state[i] * state[j];
			stateIndex++;
		}
		return secondMoments;
	}

	@Override
	protected void sampleSubnetworkState(double t, double[] x,
			SubnetworkInformation subnetworkInfo) {
		// TODO Auto-generated method stub
		
	}

}
