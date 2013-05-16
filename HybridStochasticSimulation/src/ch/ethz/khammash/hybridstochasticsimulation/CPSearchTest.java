package ch.ethz.khammash.hybridstochasticsimulation;

import JaCoP.core.*; 
import JaCoP.constraints.*; 
import JaCoP.search.*; 
 
public class CPSearchTest { 
 
    public int[] search () { 
        Store store = new Store();  // define FD store 
        int size = 4; 
        // define finite domain variables 
        IntVar[] v = new IntVar[size]; 
        for (int i=0; i<size; i++) 
            v[i] = new IntVar(store, "v"+i, 1, size); 
        // define constraints 
        store.impose( new XneqY(v[0], v[1]) ); 
        store.impose( new XneqY(v[0], v[2]) ); 
        store.impose( new XneqY(v[1], v[2]) ); 
        store.impose( new XneqY(v[1], v[3]) ); 
        store.impose( new XneqY(v[2], v[3]) ); 

        // search for a solution and print results 
        Search<IntVar> search = new DepthFirstSearch<IntVar>(); 
        SelectChoicePoint<IntVar> select = 
        new InputOrderSelect<IntVar>(store, v, 
                                 new IndomainMin<IntVar>()); 
        boolean result = search.labeling(store, select); 

        if ( result ) 
            System.out.println("Solution: " + v[0]+", "+v[1] +", "+ 
                                          v[2] +", "+v[3]); 
        else 
            System.out.println("*** No"); 

        int[] x = new int[size];
        for (int i=0; i<size; i++)
            x[i] = v[i].value();
        return x;
    }

}
