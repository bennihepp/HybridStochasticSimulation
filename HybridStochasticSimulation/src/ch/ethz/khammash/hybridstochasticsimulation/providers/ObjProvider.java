package ch.ethz.khammash.hybridstochasticsimulation.providers;


public interface ObjProvider<T> extends javax.inject.Provider<T> {

	@Override
	T get();

}
