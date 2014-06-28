package ch.ethz.bhepp.hybridstochasticsimulation.simulators;

import java.util.Collection;


public class NonEmptyFilter<T, E extends Collection<T>> extends MinimumSizeFilter<T, E> {

	public NonEmptyFilter() {
		super(1);
	}

}
