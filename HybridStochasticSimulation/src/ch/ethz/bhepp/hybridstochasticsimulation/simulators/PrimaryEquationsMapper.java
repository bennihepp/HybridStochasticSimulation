package ch.ethz.bhepp.hybridstochasticsimulation.simulators;

import org.apache.commons.math3.exception.DimensionMismatchException;
import org.apache.commons.math3.ode.EquationsMapper;


public class PrimaryEquationsMapper extends EquationsMapper {

	private static final long serialVersionUID = -3648083568032591276L;

	private final int dimension;

    public PrimaryEquationsMapper(final int dimension) {
    	super(0, dimension);
        this.dimension  = dimension;
    }

    /** Get the index of the first equation element in complete state arrays.
     * @return index of the first equation element in complete state arrays
     */
    @Override
    public int getFirstIndex() {
        return 0;
    }

    /** Get the dimension of the secondary state parameters.
     * @return dimension of the secondary state parameters
     */
    @Override
    public int getDimension() {
        return dimension;
    }

    /** Extract equation data from a complete state or derivative array.
     * @param complete complete state or derivative array from which
     * equation data should be retrieved
     * @param equationData placeholder where to put equation data
     * @throws DimensionMismatchException if the dimension of the equation data does not
     * match the mapper dimension
     */
    public void extractEquationData(double[] complete, double[] equationData)
        throws DimensionMismatchException {
//        if (equationData.length != dimension) {
//            throw new DimensionMismatchException(equationData.length, dimension);
//        }
//    	System.out.println("extractEquationData");
        System.arraycopy(complete, 0, equationData, 0, dimension);
    }

    /** Insert equation data into a complete state or derivative array.
     * @param equationData equation data to be inserted into the complete array
     * @param complete placeholder where to put equation data (only the
     * part corresponding to the equation will be overwritten)
     * @throws DimensionMismatchException if the dimension of the equation data does not
     * match the mapper dimension
     */
    public void insertEquationData(double[] equationData, double[] complete)
        throws DimensionMismatchException {
//        if (equationData.length != dimension) {
//            throw new DimensionMismatchException(equationData.length, dimension);
//        }
//    	System.out.println("insertEquationData");
        System.arraycopy(equationData, 0, complete, 0, dimension);
    }

}
