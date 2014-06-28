package ch.ethz.bhepp.hybridstochasticsimulation.averaging;

import java.util.Iterator;

import com.google.common.base.Predicate;
import com.google.common.collect.Iterators;

public class FilteredSubnetworkEnumerator implements SubnetworkEnumerator {

	private SubnetworkEnumerator baseEnumerator;
	private Predicate<SubnetworkDescription> filter;

	public FilteredSubnetworkEnumerator(SubnetworkEnumerator baseEnumerator, Predicate<SubnetworkDescription> filter) {
		this.baseEnumerator = baseEnumerator;
		this.filter = filter;
	}

	@Override
	public Iterator<SubnetworkDescription> iterator() {
		Iterator<SubnetworkDescription> baseIterator = baseEnumerator.iterator();
		return Iterators.filter(baseIterator, filter);
	}

}
