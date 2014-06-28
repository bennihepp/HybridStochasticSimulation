package ch.ethz.bhepp.hybridstochasticsimulation.providers;


public interface ObjProvider<T> extends javax.inject.Provider<T> {

	@Override
	T get();

}
