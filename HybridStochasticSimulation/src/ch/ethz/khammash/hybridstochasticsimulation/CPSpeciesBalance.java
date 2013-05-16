package ch.ethz.khammash.hybridstochasticsimulation;

import JaCoP.core.*;
import JaCoP.constraints.*;
import JaCoP.search.*;
import java.util.Vector;

public class CPSpeciesBalance {

    int[][] stoch_matrix;

    static int INT_MIN = -1000;
    static int INT_MAX = 1000;

    public static void main(String[] args) {
    	CPSpeciesBalance cpsb = new CPSpeciesBalance();
    	int r0 = 2;
    	int s0 = 4;
    	int[][] stoch_matrix = new int[r0][s0];
    	stoch_matrix[0][0] = -2;
    	stoch_matrix[0][1] = 1;
    	stoch_matrix[0][2] = 0;
    	stoch_matrix[0][3] = 0;
    	stoch_matrix[1][0] = -1;
    	stoch_matrix[1][1] = 0;
    	stoch_matrix[1][2] = -1;
    	stoch_matrix[1][3] = 1;
    	cpsb.setStochMatrix(stoch_matrix);

    	int[] min_alpha = {1, 0, 0, 0};
    	int[] max_alpha = {1, 10, 10, 10};
    	int[] min_beta = {-10, -10};
    	int[] max_beta = {10, 10};
    	double gamma = 0.0;

    	int[] x = cpsb.search_solution(min_alpha, max_alpha, min_beta, max_beta, gamma);
    	for (int i=0; i < x.length; i++)
    		System.out.println("x["+i+"] = "+x[i]);
    }

    public CPSpeciesBalance() {
        this.stoch_matrix = null;
    }

    public int getNumberOfSpecies() {
    	return this.stoch_matrix[0].length;
    }
    public int getNumberOfReactions() {
    	return this.stoch_matrix.length;
    }

    public void setStochMatrixElement(int r, int s, int value) {
		this.stoch_matrix[r][s] = value;
    }

    public void setStochMatrix(int[][] stoch_matrix) {
    	int r0 = stoch_matrix.length;
    	int s0 = stoch_matrix[0].length;
    	this.stoch_matrix = new int[r0][s0];
    	for (int r=0; r < r0; r++)
    		for (int s=0; s < s0; s++)
    			this.stoch_matrix[r][s] = stoch_matrix[r][s];
    }

    public int[] search_solution(
    		int[] min_alpha, int[] max_alpha,
    		int[] min_beta, int[] max_beta,
    		double gamma) {
        // Define FD store
        Store store = new Store();
        int s0 = min_alpha.length;
        int r0 = min_beta.length;
        // Define finite domain variables
        IntVar[] alpha = new IntVar[s0];
        for (int i=0; i < min_alpha.length; i++)
            alpha[i] = new IntVar(store, "alpha"+i,
                min_alpha[i], max_alpha[i]);
        IntVar[] beta = new IntVar[r0];
        for (int i=0; i < min_beta.length; i++)
            beta[i] = new IntVar(store, "beta"+i,
                min_beta[i], max_beta[i]);
        IntVar[] rho = build_rho(store, alpha, beta, gamma);

        Vector<PrimitiveConstraint[]> constraints = build_species_balance_constraints(
        		store, alpha, beta, rho, gamma);

        // Try balance equation constraints first
        for (int s=0; s < s0; s++) {
        	Constraint bs_constraint = constraints.get(s)[0];
        	Constraint ts_constraint = constraints.get(s)[1];
        	if (bs_constraint != null)
        		store.impose(bs_constraint);
        	else if (ts_constraint != null)
        		store.impose(ts_constraint);
        }

        IntVar[] v_array = new IntVar[s0 + r0 + 1];
        int k = 0;
        for (int i=0; i < min_alpha.length; i++)
            v_array[k++] = alpha[i];
        for (int i=0; i < min_beta.length; i++)
            v_array[k++] = beta[i];

        // search for a solution and print results
        Search<IntVar> search = new DepthFirstSearch<IntVar>();
        SelectChoicePoint<IntVar> select =
        		new InputOrderSelect<IntVar>(store, v_array,
        				new IndomainMiddle<IntVar>());
        boolean result = search.labeling(store, select);

        if (result) {
            System.out.println("Solution found");
            int[] x = new int[s0 + r0];
            for (int i=0; i < s0 + r0; i++)
                x[i] = v_array[i].value();
            return x;
        }

        // Now try balance equation and time scale constraint
        store = new Store();
        // Define finite domain variables
        alpha = new IntVar[s0];
        for (int i=0; i < min_alpha.length; i++)
            alpha[i] = new IntVar(store, "alpha"+i,
                min_alpha[i], max_alpha[i]);
        beta = new IntVar[r0];
        for (int i=0; i < min_beta.length; i++)
            beta[i] = new IntVar(store, "beta"+i,
                min_beta[i], max_beta[i]);
        rho = build_rho(store, alpha, beta, gamma);

        constraints = build_species_balance_constraints(
        		store, alpha, beta, rho, gamma);

        for (int s=0; s < s0; s++) {
        	PrimitiveConstraint bs_constraint = constraints.get(s)[0];
        	PrimitiveConstraint ts_constraint = constraints.get(s)[1];

            Constraint c = null;
            if (bs_constraint != null && ts_constraint != null) {
                c = new Or(bs_constraint, ts_constraint);
            } else
                c = ts_constraint;

            if (c != null)
                store.impose(c);
        }

        v_array = new IntVar[s0 + r0 + 1];
        k = 0;
        for (int i=0; i < min_alpha.length; i++)
            v_array[k++] = alpha[i];
        for (int i=0; i < min_beta.length; i++)
            v_array[k++] = beta[i];

        // search for a solution and print results
        search = new DepthFirstSearch<IntVar>();
        select =
        		new InputOrderSelect<IntVar>(store, v_array,
        				new IndomainMiddle<IntVar>());
        result = search.labeling(store, select);

        if (result)
            System.out.println("Solution found");
        else
            System.out.println("*** No solution found!");

        int[] x = new int[s0 + r0];
        for (int i=0; i < s0 + r0; i++)
            x[i] = v_array[i].value();
        return x;
    }

    IntVar[] build_rho(Store store, IntVar[] alpha, IntVar[] beta, double gamma) {
    	int s0 = alpha.length;
    	int r0 = beta.length;

        // Define rhos
        Vector<IntVar> rho_vec = new Vector<IntVar>();
        for (int r=0; r < r0; r++) {
            Vector<IntVar> vlist = new Vector<IntVar>();
            vlist.add(beta[r]);
        	for (int s=0; s < s0; s++) {
                if (this.stoch_matrix[r][s] < 0) {
                    IntVar p = new IntVar(store, INT_MIN, INT_MAX);
                    store.impose(new XmulCeqZ(
                             alpha[s], -this.stoch_matrix[r][s], p));
                    vlist.add(p);
                }
        	}
    		IntVar q = new IntVar(store, "rho"+r, INT_MIN, INT_MAX);
            store.impose(new Sum(convertVectorToArray(vlist), q));
            rho_vec.add(q);
        }
        IntVar[] rho = convertVectorToArray(rho_vec);
        return rho;
    }

    Vector<PrimitiveConstraint[]> build_species_balance_constraints(
    		Store store, IntVar[] alpha, IntVar[] beta, IntVar[] rho, double gamma) {
    	int s0 = alpha.length;
    	int r0 = beta.length;

    	Vector<PrimitiveConstraint[]> constraints = new Vector<PrimitiveConstraint[]>();

        // Define constraints
        for (int s=0; s < s0; s++) {
        	PrimitiveConstraint bs_constraint = null;
        	PrimitiveConstraint ts_constraint = null;

            // Build variables reflecting the max inside scaling
        	Vector<IntVar> producing_rho = new Vector<IntVar>();
        	for (int r=0; r < r0; r++) {
        		if (this.stoch_matrix[r][s] > 0)
        			producing_rho.add(rho[r]);
        	}
        	Vector<IntVar> consuming_rho = new Vector<IntVar>();
        	for (int r=0; r < r0; r++) {
        		if (this.stoch_matrix[r][s] < 0)
        			consuming_rho.add(rho[r]);
        	}

        	IntVar producing_rho_max = null;
        	IntVar consuming_rho_max = null;
        	if (producing_rho.size() > 0) {
	            producing_rho_max = new IntVar(store, INT_MIN, INT_MAX);
	            store.impose(new Max(convertVectorToArray(producing_rho), producing_rho_max));
        	}
        	if (consuming_rho.size() > 0) {
	            consuming_rho_max = new IntVar(store, INT_MIN, INT_MAX);
	            store.impose(new Max(convertVectorToArray(consuming_rho), consuming_rho_max));
        	}

            if (producing_rho_max != null && consuming_rho_max != null)
                bs_constraint = new XeqY(producing_rho_max, consuming_rho_max);

        	IntVar total_rho_max = null;
        	if (producing_rho_max != null && consuming_rho_max != null) {
	            total_rho_max = new IntVar(store, INT_MIN, INT_MAX);
            	IntVar[] rho_max_arr = {producing_rho_max, consuming_rho_max};
	            store.impose(new Max(rho_max_arr, total_rho_max));
        	} else if (producing_rho_max != null) {
        		total_rho_max = producing_rho_max;
        	} else if (consuming_rho_max != null) {
        		total_rho_max = consuming_rho_max;
        	}

            // check time-scale constraint
            if (total_rho_max != null)
            {
                IntVar total_rho_max_neg = new IntVar(store, INT_MIN, INT_MAX);
                store.impose(new XmulCeqZ(total_rho_max, -1, total_rho_max_neg));
                IntVar time_scale = new IntVar(store, INT_MIN, INT_MAX);
                store.impose(new XplusYeqZ(alpha[s], total_rho_max_neg, time_scale));
                int gamma_int = gamma >= 0 ? (int)gamma : (int)(gamma-1);
                ts_constraint = new XgteqC(time_scale, gamma_int);
            }

            PrimitiveConstraint[] cs = new PrimitiveConstraint[2];
            cs[0] = bs_constraint;
            cs[1] = ts_constraint;
            constraints.add(cs);
        }

        return constraints;
    }

    IntVar[] convertVectorToArray(Vector<IntVar> vec) {
        IntVar[] arr = new IntVar[vec.size()];
        for (int i=0; i < vec.size(); i++)
            arr[i] = vec.get(i);
        return arr;
    }

}
