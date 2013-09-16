package ch.ethz.khammash.hybridstochasticsimulation.averaging;

import java.util.Arrays;
import java.util.Iterator;
import java.util.List;
import java.util.Stack;

import org.apache.commons.math3.util.FastMath;
import org.ejml.alg.dense.mult.VectorVectorMult;
import org.ejml.data.DenseMatrix64F;

import com.google.common.base.Supplier;

import ch.ethz.khammash.hybridstochasticsimulation.graphs.SpeciesVertex;

public class SpeciesConservationRelationStateEnumerator implements Iterable<double[]> {

	private enum StateConsitence {
		UNDERFLOW, CONSISTENT, OVERFLOW,
	}

	private final int relationConstant;
	private final int[] maxState;
	private final int[] lc;
	private final Supplier<Integer> numberOfStatesSupplier;

	public SpeciesConservationRelationStateEnumerator(double[] x, SpeciesConservationRelation relation) {
		List<SpeciesVertex> speciesList = relation.getConservedSpeciesList();
		DenseMatrix64F lcVector = relation.getLinearCombination();
		DenseMatrix64F xVector = new DenseMatrix64F(lcVector.getNumRows(), 1);
		for (int i=0; i < xVector.getNumRows(); i++)
			xVector.set(i, x[speciesList.get(i).getSpecies()]);
		// TODO: Is it ok to round here?
		relationConstant = (int)FastMath.round(VectorVectorMult.innerProd(lcVector, xVector));
		lc = new int[lcVector.getNumRows()];
		for (int i=0; i < lcVector.getNumRows(); i++)
			lc[i] = (int)FastMath.round(lcVector.get(i));
		maxState = new int[speciesList.size()];
		for (int i=0; i < maxState.length; i++)
			maxState[i] = relationConstant / lc[i];
		numberOfStatesSupplier = new Supplier<Integer>() {

			@Override
			public Integer get() {
				int numOfStates = 0;
				Iterator<double[]> it = iterator();
				while (it.hasNext())
					numOfStates++;
				return numOfStates;
			}

		};
	}

	/**
	 * Computes the number of possible states.
	 *
	 * @return the number of possible states
	 */
	// This method is final because it strictly depends on the implementation of the iterator() method
	// TODO: This implementation is inefficient (simply looping through all possible states)
	final public int getNumberOfStates() {
		return numberOfStatesSupplier.get();
	}

	/**
	 * Computes the enumeration index for the provided state.
	 *
	 * @param a valid (sum_i (state_i) = relationConstant) state (the validity is not checked!)
	 * @return the enumeration index of the state
	 */
	// This method is final because it strictly depends on the implementation of the iterator() method
	// TODO: This implementation is inefficient (simply looping through all possible states)
	final public int getStateIndex(double[] state) {
//		checkArgument(state.length == maxState.length,
//				"Expected state.length == maxState.length, got %d != %d", state.length, maxState.length);
////		int stateIndex = maxState[state.length] / lc[state.length];
//		int stateIndex = 0;
//		int cumulatedNumOfStates = 1;
//		for (int j=state.length - 2; j >= 0; j--) {
//			stateIndex += cumulatedNumOfStates * state[j];
//			cumulatedNumOfStates *= maxState[j] / lc[j];
//		}
//		return stateIndex;
		int stateIndex = 0;
		for (double[] itState : this) {
			if (Arrays.equals(state, itState))
				return stateIndex;
			stateIndex++;
		}
		throw new IllegalArgumentException("The passed state is not valid");
	}

	// This method is final because the method getStateIndex() depends depends on the iteration order
	@Override
	final public Iterator<double[]> iterator() {
		return new Iterator<double[]>() {

			private double[] state;
			private Stack<Integer> stack;
			{
				state = new double[maxState.length];
				stack = new Stack<>();
				state[state.length - 2] = -1;
				stack.push(state.length - 2);
			}

			@Override
			public boolean hasNext() {
				if (stack.isEmpty())
					return false;
				int index = stack.peek();
				state[index]++;
				boolean result = false;
				switch (validateState()) {
				case CONSISTENT:
					result = true;
					break;
				case OVERFLOW:
					state[index] = 0;
					if (index < state.length - 2)
						stack.push(index + 1);
					else
						stack.pop();
					result = hasNext();
					break;
				case UNDERFLOW:
					result = hasNext();
					break;
				}
				return result;
			}

			private StateConsitence validateState() {
				int sum = 0;
				for (int i=0; i < state.length - 1; i++)
					sum += state[i] * lc[i];
				int residual = relationConstant - sum;
				if (residual < 0)
					return StateConsitence.OVERFLOW;
				if (residual % lc[state.length - 1] == 0) {
					state[state.length - 1] = residual / lc[state.length - 1];
					return StateConsitence.CONSISTENT;
				} else
					return StateConsitence.UNDERFLOW;
			}

			@Override
			public double[] next() {
				return state.clone();
			}

			@Override
			public void remove() {
				throw new UnsupportedOperationException("This iterator doesn't support remove operations");
			}

		};
	}

}
