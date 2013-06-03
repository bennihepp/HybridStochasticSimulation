package ch.ethz.khammash.hybridstochasticsimulation.plotting;

import org.rosuda.JRI.Rengine;

public class JRItest {

	public static void main (String[] args) {
		// new R-engine
		if (!Rengine.versionCheck())
		    throw new UnsupportedClassVersionError("JRI version mismatch");
		Rengine re=new Rengine (new String [] {"--vanilla"}, false, null);
		if (!re.waitForR()) {
			System.out.println ("Cannot load R");
			return;
		}

		// print a random number from uniform distribution
		System.out.println (re.eval ("runif(1)").asDouble ());

		// done...
		re.end();
	}

}
