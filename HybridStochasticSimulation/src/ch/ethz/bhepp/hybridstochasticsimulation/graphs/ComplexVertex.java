package ch.ethz.bhepp.hybridstochasticsimulation.graphs;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class ComplexVertex {

	private int[] stochiometries;

	public ComplexVertex(int[] stochiometries) {
		this.stochiometries = stochiometries.clone();
	}

	public List<Integer> getStochiometryList() {
		ArrayList<Integer> result = new ArrayList<Integer>(stochiometries.length);
		for (int s=0; s < stochiometries.length; s++)
			result.add(stochiometries[s]);
		return result;
	}

	public int[] getStochiometries() {
		return stochiometries.clone();
	}

	@Override
	public boolean equals(Object other) {
		if (other == null)
			return false;
		if (other == this)
			return true;
		if (!other.getClass().equals(getClass()))
			return false;
		ComplexVertex otherVertex = (ComplexVertex)other;
		return Arrays.equals(this.stochiometries, otherVertex.stochiometries);
	}

	@Override
	public int hashCode() {
		return Arrays.hashCode(stochiometries);
	}

	@Override
	public String toString() {
		StringBuilder sb = new StringBuilder();
		sb.append("{");
		for (int s=0; s < stochiometries.length; s++) {
			sb.append(stochiometries[s]);
			sb.append(", ");
		}
		if (stochiometries.length > 0)
			sb.delete(sb.length() - 2, sb.length());
		sb.append("}");
		return sb.toString();
	}

}
