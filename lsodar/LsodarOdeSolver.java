class LsodarOdeSolver<T extends PDMPModel> {

    static {
        // Load the shared library liblsodarjni
        System.loadLibrary("liblsodarjni");
    }

    public LsodarOdeSolver() {
    };

    private native void jni_initialize(LsodarOde ode);

    public void initialize(LsodarOde ode, FiniteTrajectoryRecorder<T> tr) {
        jni_initialize(ode, tr);
    }

    private native void jni_cleanup();

    public void cleanup() {
        jni_cleanup();
    }

    private native double integrate(double t0, double[] x0, double t1);

    public double integrate(double t0, double[] x0, double t1) {
        return jni_integrate(t0, x0, t1);
    }

    //private static native Object callbackTest2(double arr[]);

    //private static native int callbackTest(double arr[]);

    //public static double callback(double v) {
    //    return v + 1;
    //}

    /*public static void main(String args[]) {
        double[] q = new double[5];
        q[0] = 9;
        q[1] = 7;
        q[2] = 5;
        q[3] = 3;
        q[4] = 1;
        System.out.println(callbackTest(q));
        for (int i=0; i < q.length; i++)
            System.out.println("q["+i+"] = " + q[i]);
    };*/

};
