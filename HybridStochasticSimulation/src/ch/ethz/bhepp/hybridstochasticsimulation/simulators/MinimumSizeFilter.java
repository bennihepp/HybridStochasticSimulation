package ch.ethz.bhepp.hybridstochasticsimulation.simulators;

import java.util.Collection;

import com.google.common.base.Predicate;


public class MinimumSizeFilter<T, E extends Collection<T>> implements Predicate<E> {

	private final int minSize;

	public MinimumSizeFilter(final int minSize) {
		this.minSize = minSize;
	}

	@Override
	public boolean apply(E input) {
		return input.size() >= minSize;
	}

}
