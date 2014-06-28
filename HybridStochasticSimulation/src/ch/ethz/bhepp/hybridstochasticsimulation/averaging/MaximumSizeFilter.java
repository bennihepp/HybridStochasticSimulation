package ch.ethz.bhepp.hybridstochasticsimulation.averaging;

import java.util.Collection;

import com.google.common.base.Predicate;

public class MaximumSizeFilter<T, E extends Collection<T>> implements Predicate<E> {

	private final int maxSize;

	public MaximumSizeFilter(final int maxSize) {
		this.maxSize = maxSize;
	}

	@Override
	public boolean apply(E species) {
		return species.size() <= maxSize;
	}

}
