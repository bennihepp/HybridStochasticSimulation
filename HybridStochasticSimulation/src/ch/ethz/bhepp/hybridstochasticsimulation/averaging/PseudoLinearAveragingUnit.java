package ch.ethz.bhepp.hybridstochasticsimulation.averaging;

import java.util.Arrays;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.DecompositionSolver;
import org.apache.commons.math3.linear.LUDecomposition;
import org.apache.commons.math3.linear.QRDecomposition;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
import org.apache.commons.math3.linear.SingularValueDecomposition;
import org.apache.commons.math3.util.FastMath;
import org.ejml.data.DenseMatrix64F;
import org.ejml.ops.CommonOps;

import ch.ethz.bhepp.hybridstochasticsimulation.math.BroydenRootSolver;
import ch.ethz.bhepp.hybridstochasticsimulation.math.MultivariateFunction;
import ch.ethz.bhepp.hybridstochasticsimulation.models.AdaptiveMSHRNModel;
import ch.ethz.bhepp.hybridstochasticsimulation.models.UnaryBinaryDeterministicModel;
import ch.ethz.bhepp.hybridstochasticsimulation.models.UnaryBinaryModelUtils;
import ch.ethz.bhepp.hybridstochasticsimulation.networks.MassActionReactionNetwork;
import ch.ethz.bhepp.hybridstochasticsimulation.networks.MSHybridReactionNetwork.ReactionType;
import ch.ethz.bhepp.hybridstochasticsimulation.networks.MSHybridReactionNetwork.SpeciesType;

import com.google.common.base.Predicate;

public class PseudoLinearAveragingUnit extends AbstractAveragingUnit {

//	private static final Logger logger = LoggerFactory.getLogger(PseudoLinearAveragingUnit.class);

	private List<SubnetworkDescription> pseudoLinearSubnetworks = null;
//	private List<SubnetworkDescription> averagingCandidates;
//	private boolean averagingInvalid;
	private boolean _stopIfAveragingBecomesInvalid = true;
	private boolean _warnIfAveragingBecomesInvalid = true;
	private boolean _performPseudoLinearAveragingOnlyOnce = true;
	private Map<SubnetworkDescription, Boolean> hasContinuousSurroundingReactionMap = null;
	private Map<SubnetworkDescription, double[]> precomputedSolutionMap = null;
	private Map<SubnetworkDescription, RealMatrix> linearEquationSystemMatrixAMap = null;
	private Map<SubnetworkDescription, RealVector> linearEquationSystemVectorBMap = null;
	private double[] _currentFirstMoments = null;

	public static PseudoLinearAveragingUnit createCopy(PseudoLinearAveragingUnit provider) {
		PseudoLinearAveragingUnit copy = new PseudoLinearAveragingUnit();
		copy.copyFrom(provider);
		copy.pseudoLinearSubnetworks = provider.pseudoLinearSubnetworks;
		copy._stopIfAveragingBecomesInvalid = provider._stopIfAveragingBecomesInvalid;
		copy._warnIfAveragingBecomesInvalid = provider._warnIfAveragingBecomesInvalid;
		copy._performPseudoLinearAveragingOnlyOnce = provider._performPseudoLinearAveragingOnlyOnce;
		return copy;
	}

	public PseudoLinearAveragingUnit(MassActionReactionNetwork network, Set<Integer> importantSpecies) {
		super(network, importantSpecies);
	}

	protected PseudoLinearAveragingUnit() {
		super();
	}

	public void stopIfAveragingBecomesInvalid(boolean stop) {
		_stopIfAveragingBecomesInvalid = stop;
	}

	public void warnIfAveragingBecomesInvalid(boolean warn) {
		_warnIfAveragingBecomesInvalid = warn;
	}

	public void performPseudoLinearAveragingOnlyOnce(boolean onlyOnce) {
		_performPseudoLinearAveragingOnlyOnce = onlyOnce;
	}

//	private List<SubnetworkDescription> findPseudoLinearSubnetworks() {
//		List<SubnetworkDescription> pseudoLinearSubnetworks = new ArrayList<>();
//		for (SubnetworkDescription subnetwork : enumerateSubnetworks()) {
//			if (isPseudoLinear(subnetwork))
//				pseudoLinearSubnetworks.add(subnetwork);
//		}
//		return pseudoLinearSubnetworks;
//	}

	private boolean isPseudoLinear(SubnetworkDescription subnetwork) {
		Set<Integer> involvedReactions = findInvolvedReactions(subnetwork.getSubnetworkSpecies());
		// Check if all involved reactions are pseudo-linear with respect to the subnetwork
		for (int r : involvedReactions) {
			int[] choiceIndices = network.getReactantIndices(r);
			int c = 0;
			for (int s : choiceIndices) {
				if (subnetwork.getSubnetworkSpecies().contains(s))
					c++;
			}
			if (c >= 2)
				return false;
		}
//		SubnetworkInformation subnetworkInfo = computeSubnetworkInformation(subnetwork);
		MassActionReactionNetwork subReactionNetwork = createSubReactionNetwork(network, subnetwork.getSubnetworkSpecies(), subnetwork.getSubnetworkReactions());
		RealMatrix a = UnaryBinaryModelUtils.computeUnaryCoefficientMatrix(subReactionNetwork);
		RealMatrix c = UnaryBinaryModelUtils.computeAugmentedCoefficientMatrix(subReactionNetwork);
		SingularValueDecomposition svdA = new SingularValueDecomposition(a);
		int rankA = svdA.getRank();
		SingularValueDecomposition svdC = new SingularValueDecomposition(c);
		int rankC = svdC.getRank();
		// Make sure that there is a stationary solution
		if (rankA != rankC)
			return false;
		return true;
	}

	private Set<Integer> findInvolvedReactions(Collection<Integer> subnetworkSpecies) {
		HashSet<Integer> involvedReactions = new HashSet<Integer>();
		for (int s : subnetworkSpecies) {
			involvedReactions.addAll(network.getInvolvingReactions(s));
		}
		return involvedReactions;
	}

	@Override
	public void reset() {
//		averagingCandidates = null;
//		averagingInvalid = false;
		pseudoLinearSubnetworks = null;
		hasContinuousSurroundingReactionMap = new HashMap<>();
		precomputedSolutionMap = new HashMap<>();
		linearEquationSystemMatrixAMap = new HashMap<>();
		linearEquationSystemVectorBMap = new HashMap<>();
		super.reset();
	}

//	@Override
//	public List<SubnetworkDescription> findAveragingCandidates(double t, double[] x, Predicate<SubnetworkDescription> filter) {
//		if (pseudoLinearSubnetworks == null)
//			pseudoLinearSubnetworks = findPseudoLinearSubnetworks();
//		if (!_performPseudoLinearAveragingOnlyOnce)
//			averagingCandidates = null;
//		if (averagingCandidates == null) {
//			averagingCandidates = new ArrayList<>();
//			for (SubnetworkDescription subnetwork : pseudoLinearSubnetworks) {
//				if (filter.apply(subnetwork))
//					averagingCandidates.add(subnetwork);
////				if (checkAveragingConditions(subnetwork, x, reactionTimescales))
////					averagingCandidates.add(subnetwork);
//			}
//		} else {
//			// Check whether the conditions for averaging of the pseudo linear subnetworks are still valid
//			for (SubnetworkDescription subnetwork : averagingCandidates) {
//				boolean satisfied = filter.apply(subnetwork);
////				boolean satisfied = checkAveragingConditions(subnetwork, x, reactionTimescales);
//				if (satisfied && averagingInvalid) {
//					averagingInvalid = false;
//					if (_warnIfAveragingBecomesInvalid)
//						logger.warn(getRevalidatedWarningMessage(t));
//				}
//				if (!satisfied && !averagingInvalid) {
//					averagingInvalid = true;
//					if (_stopIfAveragingBecomesInvalid)
//						throw new AveragingException(getInvalidWarningMessage(t));
//					if (_warnIfAveragingBecomesInvalid)
//						logger.warn(getInvalidWarningMessage(t));
//				}
//			}
//		}
//		return averagingCandidates;
//	}

//	private String getInvalidWarningMessage(double t) {
//		return "Averaging of pseudo linear subnetworks isn't valid anymore at t=" + t;
//	}

//	private String getRevalidatedWarningMessage(double t) {
//		return "Averaging of pseudo linear subnetworks switched to being valid again at t=" + t;
//	}

//	@Override
//	protected void computeAverageStationaryStateOfSubnetworks(double t, double[] x, List<SubnetworkDescription> subnetworksToAverage) {
//		// TODO: Check if this is ok.
//		// As the stationary state of a pseudo-linear subnetwork is the same as the stationary solution of the corresponding
//		// deterministic system, we don't compute the stationary state and let the PDMP solver evolve the correct solution.
//	}

	@Override
	public void updateAveraging(AdaptiveMSHRNModel model, double t, double[] x, SubnetworkDescription subnetworkDescr) {
		SubnetworkInformation subnetworkInfo = super.createSubnetworkInformation(subnetworkDescr);
		computeSubnetworkInformation(subnetworkInfo);
		MassActionReactionNetwork subReactionNetwork = subnetworkInfo.getReactionNetwork();
		// TODO: this is duplicate
		RealMatrix a = UnaryBinaryModelUtils.computeUnaryCoefficientMatrix(subnetworkInfo.getReactionNetwork());
		SingularValueDecomposition svdA = new SingularValueDecomposition(a);
		int rankA = svdA.getRank();
		if (rankA == subReactionNetwork.getNumberOfSpecies()) {
			// There exists exactly one solution, so we can precompute it
			RealMatrix matrixA = UnaryBinaryModelUtils.computeUnaryCoefficientMatrix(subReactionNetwork);
			RealVector vectorB = UnaryBinaryModelUtils.computeConstitutiveCoefficientVector(subReactionNetwork).mapMultiplyToSelf(-1);

			DecompositionSolver solver = new LUDecomposition(matrixA).getSolver();
			RealVector stationaryStateVector = solver.solve(vectorB);
			double[] firstMoments = stationaryStateVector.toArray();
			precomputedSolutionMap.put(subnetworkInfo, firstMoments);

//			DecompositionSolver solver = new LUDecomposition(matrixA).getSolver();
//			RealVector stationaryStateVector = solver.solve(vectorB);
//			precomputedSolutionMap.put(subnetworkInfo, stationaryStateVector.toArray());
		} else {
			// If there are infinitely many solutions a conservation relation has to be incorporated
			Collection<SpeciesConservationRelation> conservedSpeciesRelations = subnetworkInfo.getConservedSpeciesRelations();
			RealMatrix matrixA = UnaryBinaryModelUtils.computeUnaryCoefficientMatrix(subReactionNetwork);
			RealMatrix extMatrixA = new Array2DRowRealMatrix(
					subReactionNetwork.getNumberOfSpecies() + conservedSpeciesRelations.size(), subReactionNetwork.getNumberOfSpecies());
			for (int r=0; r < subReactionNetwork.getNumberOfSpecies(); r++)
				extMatrixA.setRow(r, matrixA.getRow(r));
			int i = subReactionNetwork.getNumberOfSpecies();
			for (SpeciesConservationRelation relation : conservedSpeciesRelations) {
				for (int j=0; j < relation.getConservedSpeciesList().size(); j++) {
					int s = relation.getConservedSpeciesList().get(j);
					extMatrixA.setEntry(i, s, relation.getLinearCombination().get(j));
				}
				i++;
			}
			RealVector vectorB = UnaryBinaryModelUtils.computeConstitutiveCoefficientVector(subReactionNetwork);
			RealVector extVectorB = new ArrayRealVector(extMatrixA.getRowDimension());
			for (int r=0; r < vectorB.getDimension(); r++)
				extVectorB.setEntry(r, vectorB.getEntry(r));
			linearEquationSystemMatrixAMap.put(subnetworkInfo, extMatrixA);
			linearEquationSystemVectorBMap.put(subnetworkInfo, extVectorB);
		}

		hasContinuousSurroundingReactionMap.put(subnetworkDescr, false);
		for (int r : subnetworkDescr.getSurroundingReactions())
			if (model.getNetwork().getReactionType(r) == ReactionType.CONTINUOUS) {
				hasContinuousSurroundingReactionMap.put(subnetworkDescr, true);
				break;
			}

		if (hasContinuousSurroundingReactionMap.get(subnetworkDescr)) {
			computeAndStoreSubnetworkInformation(subnetworkDescr);
	
			for (int r : subnetworkDescr.getSubnetworkReactions()) {
				model.getNetwork().overrideReactionType(r, ReactionType.CONTINUOUS);
			}
			for (int s : subnetworkDescr.getSubnetworkSpecies()) {
				model.getNetwork().overrideSpeciesType(s, SpeciesType.CONTINUOUS);
			}
			model.update();

			updateSubnetworkState(model, t, x, subnetworkDescr);
		} else
			super.updateAveraging(model, t, x, subnetworkDescr);
	}

	@Override
	public void updateSubnetworkState(AdaptiveMSHRNModel model, double t, double[] x, SubnetworkDescription subnetwork) {
		if (hasContinuousSurroundingReactionMap.get(subnetwork))
			return;

		SubnetworkInformation subnetworkInfo = subnetworkInformationMap.get(subnetwork);

		_currentFirstMoments = precomputedSolutionMap.get(subnetworkInfo);
		if (_currentFirstMoments == null)
			_currentFirstMoments = computeStationaryFirstMoments(t, x, subnetworkInfo);
		super.updateSubnetworkState(model, t, x, subnetworkInfo);
		_currentFirstMoments = null;
	}

	@Override
	protected double[] computeStationaryFirstMoments(double t, double[] x, SubnetworkInformation subnetworkInfo) {
		if (_currentFirstMoments != null)
			return _currentFirstMoments;

		Map<Integer, Integer> indexMap = subnetworkInfo.getIndexMap();

//		RealVector vectorX = new ArrayRealVector(x, false);
		RealMatrix matrixA = linearEquationSystemMatrixAMap.get(subnetworkInfo);
		RealVector vectorB = linearEquationSystemVectorBMap.get(subnetworkInfo);

		for (int r=subnetworkInfo.getSubnetworkSpecies().size(); r < matrixA.getRowDimension(); r++) {
			double v = 0.0;
			for (int i=0; i < matrixA.getColumnDimension(); i++)
				v += matrixA.getEntry(r, i) * x[indexMap.get(i)];
//			double v = matrixA.getRowVector(r).dotProduct(vectorX);
			vectorB.setEntry(r, v);
		}

		DecompositionSolver solver = new QRDecomposition(matrixA).getSolver();
		RealVector stationaryStateVector = solver.solve(vectorB);
		double[] firstMoments = stationaryStateVector.toArray();
		return firstMoments;

//		double[] firstMoments = precomputedSolutionMap.get(subnetworkInfo);
//		if (firstMoments != null)
//			return firstMoments;
//
//		firstMoments = subnetworkInfo.computeStationarySolution(t, x);
//		UnaryBinaryDeterministicModel subnetworkModel = subnetworkInfo.getDeterministicModel();
//		Map<Integer, Integer> indexMap = subnetworkInfo.getIndexMap();
//
//		double[] initialX = new double[subnetworkModel.getNumberOfSpecies()];
//		for (int s : subnetworkInfo.getSubnetworkSpecies()) {
//			initialX[indexMap.get(s)] = x[s];
//		}
//
//		firstMoments = UnaryBinaryModelUtils.computeStationarySolution(subnetworkModel, t, initialX);
//
//		for (int s=0; s < subnetworkModel.getNumberOfSpecies(); s++) {
//			if (firstMoments[s] < 0.0)
//				firstMoments[s] = 0.0;
//		}
//
//		return firstMoments;
	}

	@Override
	protected double[][] computeStationarySecondMoments(double t, final double[] x, SubnetworkInformation subnetworkInfo) {
		final double[] firstMoments = computeStationaryFirstMoments(t, x, subnetworkInfo);
		// Compute stationary solution of moment equations
		final MassActionReactionNetwork subReactionNetwork = subnetworkInfo.getReactionNetwork();
		final UnaryBinaryDeterministicModel subnetworkModel = subnetworkInfo.getDeterministicModel();
		Map<Integer, Integer> indexMap = subnetworkInfo.getIndexMap();

		DenseMatrix64F initialXX = new DenseMatrix64F(subReactionNetwork.getNumberOfSpecies(), subReactionNetwork.getNumberOfSpecies());
		for (int i : subnetworkInfo.getSubnetworkSpecies())
			for (int j : subnetworkInfo.getSubnetworkSpecies()) {
				int row = indexMap.get(i);
				int col = indexMap.get(j);
				double v = firstMoments[i] * firstMoments[j];
				initialXX.set(row, col, v);
			}

//		double[][] secondMoments = UnaryBinaryModelUtils.computeStationarySecondMoments(subnetworkModel, t, initialX);

		MultivariateFunction rootFunction = new MultivariateFunction() {

			private MassActionReactionNetwork network;
			private UnaryBinaryDeterministicModel model;
			private DenseMatrix64F S;
			private DenseMatrix64F W;
			private DenseMatrix64F w0;
			private DenseMatrix64F Q;
			private DenseMatrix64F P;
			private DenseMatrix64F R;
			private DenseMatrix64F D;
			private DenseMatrix64F T;
			private DenseMatrix64F X;
			private DenseMatrix64F XX;
			private DenseMatrix64F dXX;
			{
				network = subReactionNetwork;
				model = subnetworkModel;
				S = computeStoichiometryMatrix();
				W = computeLinearPropensityMatrix();
				w0 = computeConstantPropensityVector();
				Q = new DenseMatrix64F(network.getNumberOfSpecies(), network.getNumberOfSpecies());
				P = new DenseMatrix64F(network.getNumberOfSpecies(), 1);
				R = new DenseMatrix64F(network.getNumberOfReactions(), 1);
				D = new DenseMatrix64F(network.getNumberOfReactions(), network.getNumberOfReactions());
				T = new DenseMatrix64F(network.getNumberOfSpecies(), network.getNumberOfReactions());
				X = new DenseMatrix64F(network.getNumberOfSpecies(), 1, true, firstMoments);
				XX = new DenseMatrix64F(network.getNumberOfSpecies(), network.getNumberOfSpecies());
				dXX = new DenseMatrix64F(network.getNumberOfSpecies(), network.getNumberOfSpecies());
			}

			private DenseMatrix64F computeStoichiometryMatrix() {
				DenseMatrix64F S = new DenseMatrix64F(network.getNumberOfSpecies(), network.getNumberOfReactions());
				for (int s=0; s < S.numRows; s++)
					for (int r=0; r < S.numCols; r++)
						S.set(s, r, network.getStoichiometry(s, r));
				return S;
			}

			private DenseMatrix64F computeLinearPropensityMatrix() {
				DenseMatrix64F W = new DenseMatrix64F(network.getNumberOfReactions(), network.getNumberOfSpecies());
				for (int r=0; r < network.getNumberOfReactions(); r++) {
					int[] choiceIndices = network.getReactantIndices(r);
					if (choiceIndices.length == 2)
						throw new UnsupportedOperationException("The network is not linear");
					if (choiceIndices.length == 1)
						W.set(r, choiceIndices[0], network.getRateParameter(r));
				}
				return W;
			}

			private DenseMatrix64F computeConstantPropensityVector() {
				DenseMatrix64F w0 = new DenseMatrix64F(network.getNumberOfReactions(), 1);
				for (int r=0; r < network.getNumberOfReactions(); r++) {
					int[] choiceIndices = network.getReactantIndices(r);
					if (choiceIndices.length == 0)
						w0.set(r, 1, network.getRateParameter(r));
				}
				return w0;
			}

			@Override
			public int getDimension() {
				return model.getNumberOfSpecies() * model.getNumberOfSpecies();
			}

			@Override
			public void computeValue(double[] x, double[] y) {
				Arrays.fill(y, 0.0);
				dXX.setData(y);
				XX.setData(x);
				CommonOps.mult(S, W, Q);
				CommonOps.mult(S, w0, P);
				CommonOps.multAdd(Q, XX, dXX);
				CommonOps.multAddTransB(XX, Q, dXX);
				CommonOps.multAddTransB(P, X, dXX);
				CommonOps.multAddTransB(X, P, dXX);
				CommonOps.mult(W, X, R);
				CommonOps.add(R, w0, R);
				for (int r=0; r < network.getNumberOfReactions(); r++) {
					D.set(r, r, R.get(r, 0));
				}
				CommonOps.mult(S, D, T);
				CommonOps.multAddTransB(T, S, dXX);
			}

		};

		BroydenRootSolver solver = new BroydenRootSolver(rootFunction);
		double[] initialXXArray = initialXX.getData();
		double[] solution = solver.findRoot(initialXXArray);

		DenseMatrix64F XX = new DenseMatrix64F(subReactionNetwork.getNumberOfSpecies(), subReactionNetwork.getNumberOfSpecies(), true, solution);
		double[][] secondMoments = new double[subReactionNetwork.getNumberOfSpecies()][subReactionNetwork.getNumberOfSpecies()];
		for (int i=0; i < secondMoments.length; i++)
			for (int j=0; j < secondMoments.length; j++)
				secondMoments[i][j] = XX.get(i, j);
		return secondMoments;
	}

//	@Override
//	protected double[] computeStationaryFirstMoments(double t, double[] x, SubnetworkDescription subnetwork) {
//		// NOTE: We approximate the exact distribution of a reducible zero-deficiency subnetwork
//		// with a multinomial distribution. Thus the average stationary state is equal to the deterministic stationary state.
//		double[] firstMoments = new double[x.length];
//		SubnetworkInformation subnetworkInfo = subnetworkInformationMap.get(subnetwork);
//		UnaryBinaryReactionNetwork subReactionNetwork = subnetworkInfo.getReactionNetwork();
////		UnaryBinaryDeterministicModel subnetworkModel = subnetworkInfo.getModel();
//		Map<Integer, Integer> indexMap = subnetworkInfo.getIndexMap();
//		double[] initialX = new double[subReactionNetwork.getNumberOfSpecies()];
//		for (int s : subnetwork.getSubnetworkSpecies()) {
//			initialX[indexMap.get(s)] = x[s];
//		}
//
//		RealVector vectorInitialX = new ArrayRealVector(initialX, false);
//		RealMatrix matrixA = linearEquationSystemMatrixAMap.get(subnetwork);
//		RealVector vectorB = linearEquationSystemVectorBMap.get(subnetwork);
//		for (int r=subReactionNetwork.getNumberOfSpecies(); r < matrixA.getRowDimension(); r++) {
//			double v = matrixA.getRowVector(r).dotProduct(vectorInitialX);
//			vectorB.setEntry(r, v);
//		}
//		DecompositionSolver solver = new QRDecomposition(matrixA).getSolver();
//		RealVector subStationarySolution = solver.solve(vectorB);
//
////		RealMatrix a = UnaryBinaryModelUtils.computeUnaryCoefficientMatrix(subReactionNetwork);
////		RealVector b = UnaryBinaryModelUtils.computeConstitutiveCoefficientVector(subReactionNetwork).mapMultiplyToSelf(-1);
////		DecompositionSolver solver = new LUDecomposition(a).getSolver();
////		RealVector subXSteadyState = solver.solve(b);
//
////		double[] subXSteadyState = UnaryBinaryModelUtils.computeSteadyState(subnetworkModel, t, initialX);
//
//		for (int s : subnetwork.getSubnetworkSpecies()) {
////			firstMoments[s] = subXSteadyState[indexMap.get(s)];
//			firstMoments[s] = subStationarySolution.getEntry(indexMap.get(s));
//			// TODO: check some tolerance of negative values (firstMoment[s] < -tolerance
//			if (firstMoments[s] < 0.0)
//				firstMoments[s] = 0.0;
//		}
//		return firstMoments;
//	}

//	@Override
//	public double[][] computeStationarySecondMoments(double t, double[] x, SubnetworkDescription subnetwork) {
//		return null;
////		throw new UnsupportedOperationException();
//	}

	@Override
	public void sampleSubnetworkState(double t, double[] x, SubnetworkInformation subnetworkInfo) {
		// We don't know the stationary distribution of the subnetwork in general.
		// Instead of sampling from the stationary distribution we just round the copy numbers of the stationary first moments
		double[] firstMoments = computeStationaryFirstMoments(t, x, subnetworkInfo);
		for (int s : subnetworkInfo.getSubnetworkSpecies()) {
			double value = firstMoments[s];
			value = FastMath.round(value);
			x[s] = value;
		}
	}

	@Override
	public Predicate<SubnetworkDescription> getSubnetworkFilter() {
		return new Predicate<SubnetworkDescription>() {

			@Override
			public boolean apply(SubnetworkDescription subnetworkDescr) {
				return isPseudoLinear(subnetworkDescr);
			}

		};
	}

}
