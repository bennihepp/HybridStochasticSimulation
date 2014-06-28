package ch.ethz.bhepp.hybridstochasticsimulation.models;

import java.util.EventListener;


public interface ReactionRateBoundEventListener extends EventListener {

    void reactionRateBoundEventOccured(double t, double[] x);

    void reactionRateBoundEventOccured(int species, double t, double[] x);

}
