package ch.ethz.bhepp.hybridstochasticsimulation.simulators;

import com.google.common.base.Predicate;
import com.google.common.collect.Iterables;

public class FilterUtilities {

	public static <T> Predicate<T> and(final Predicate<T> filter1, final Predicate<T> filter2) {
		return new Predicate<T>() {

			@Override
			public boolean apply(T input) {
				return filter1.apply(input) && filter2.apply(input);
			}

		};
	}

	public static <T> Iterable<T> filter(Iterable<T> iterable, Predicate<T> filter) {
		return Iterables.filter(iterable, filter);
	}

}
