package ch.ethz.khammash.hybridstochasticsimulation.examples;

import org.apache.commons.math3.random.MersenneTwister;
import org.apache.commons.math3.random.RandomDataGenerator;
import org.apache.commons.math3.util.FastMath;

import ch.ethz.khammash.hybridstochasticsimulation.networks.DefaultUnaryBinaryReactionNetwork;



public class RepressedBacteriumOperonNetwork extends ExampleNetwork {

	public RepressedBacteriumOperonNetwork() {
		// See Crudu et. al 2009
		// species:
		//  R
		//  RNAP
		//  D
		//  D.R
		//  D.RNAP
		//  TrRNAP
		//  RBS
		//  Rib
		//  Rib.RBS
		//  ElRib
		//  Protein
		//  FoldedProtein
		// reactions:
		//  D + R -> D.R (k1)
		//  D.R -> D + R (km1)
		//  D + RNAP -> D.RNAP (k2)
		//  D.RNAP -> D + RNAP (km2)
		//  D.RNAP -> TrRNAP (k3)
		//  TrRNAP -> RBS + D + RNAP (k4)
		//  RBS -> - (k5)
		//  Rib + RBS -> Rib.RBS (k6)
		//  Rib.RBS -> Rib + RBS (km6)
		//  Rib.RBS -> ElRib + RBS (k7)
		//  ElRib -> Protein (k8)
		//  Protein -> FoldedProtein (k9)
		//  Protein -> - (k10)
		//  FoldedProtein -> - (k11)
		double[] x0 = {
				2500, // R
				35, // RNAP
				20, // D
				0, // D.R
				0, // D.RNAP
				0, // TrRNAP
				0, // RBS
				350, // Rib
				0, // Rib.RBS
				0, // ElRib
				0, // Protein
				0, // FoldedProtein
		};
		double k1 =  1e8;
		double km1 = 1;
		double k2 =  1e8;
		double km2 = 10;
		double k3 =  0.1;
		double k4 =  0.3;
		double k5 =  0.3;
		double k6 =  1e8;
		double km6 = 2.25;
		double k7 =  0.5;
		double k8 =  0.015;
		double k9 =  FastMath.log(2) / 5400.0;
		double k10 = 1e-5;
		double k11 = 1e-5;
		double[] rateParameters = {
				k1,
				km1,
				k2,
				km2,
				k3,
				k4,
				k5,
				k6,
				km6,
				k7,
				k8,
				k9,
				k10,
				k11,
		};
		double[] plotScales = { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 };
		//  D + R -> D.R (k1)
		//  D.R -> D + R (km1)
		//  D + RNAP -> D.RNAP (k2)
		//  D.RNAP -> D + RNAP (km2)
		//  D.RNAP -> TrRNAP (k3)
		//  TrRNAP -> RBS + D + RNAP (k4)
		//  RBS -> - (k5)
		//  Rib + RBS -> Rib.RBS (k6)
		//  Rib.RBS -> Rib + RBS (km6)
		//  Rib.RBS -> ElRib + RBS (k7)
		//  ElRib -> Protein (k8)
		//  Protein -> FoldedProtein (k9)
		//  Protein -> - (k10)
		//  FoldedProtein -> - (k11)
		int[][] productionStochiometries = {
		//        R  RNAP  D  D.R D.RNAP TrRNAP RBS Rib Rib.RBS ElRib Protein FoldedProtein
				{ 0,  0,   0,  1,    0,     0,   0,  0,    0,     0,     0,        0 }, // D + R -> D.R (k1)
				{ 1,  0,   1,  0,    0,     0,   0,  0,    0,     0,     0,        0 }, // D.R -> D + R (km1)
				{ 0,  0,   0,  0,    1,     0,   0,  0,    0,     0,     0,        0 }, // D + RNAP -> D.RNAP (k2)
				{ 0,  1,   1,  0,    0,     0,   0,  0,    0,     0,     0,        0 }, // D.RNAP -> D + RNAP (km2)
				{ 0,  0,   0,  0,    0,     1,   0,  0,    0,     0,     0,        0 }, // D.RNAP -> TrRNAP (k3)
				{ 0,  1,   1,  0,    0,     0,   1,  0,    0,     0,     0,        0 }, // TrRNAP -> RBS + D + RNAP (k4)
				{ 0,  0,   0,  0,    0,     0,   0,  0,    0,     0,     0,        0 }, // RBS -> - (k5)
				{ 0,  0,   0,  0,    0,     0,   0,  0,    1,     0,     0,        0 }, // Rib + RBS -> Rib.RBS (k6)
				{ 0,  0,   0,  0,    0,     0,   1,  1,    0,     0,     0,        0 }, // Rib.RBS -> Rib + RBS (km6)
				{ 0,  0,   0,  0,    0,     0,   1,  0,    0,     1,     0,        0 }, // Rib.RBS -> ElRib + RBS (k7)
				{ 0,  0,   0,  0,    0,     0,   0,  1,    0,     0,     1,        0 }, // ElRib -> Rib + Protein (k8)
				{ 0,  0,   0,  0,    0,     0,   0,  0,    0,     0,     0,        1 }, // Protein -> FoldedProtein (k9)
				{ 0,  0,   0,  0,    0,     0,   0,  0,    0,     0,     0,        0 }, // Protein -> - (k10)
				{ 0,  0,   0,  0,    0,     0,   0,  0,    0,     0,     0,        0 }, // FoldedProtein -> - (k11)
		};
		int[][] consumptionStochiometries = {
		//        R  RNAP  D  D.R D.RNAP TrRNAP RBS Rib Rib.RBS ElRib Protein FoldedProtein
				{ 1,  0,   1,  0,    0,     0,   0,  0,    0,     0,     0,        0 }, // D + R -> D.R (k1)
				{ 0,  0,   0,  1,    0,     0,   0,  0,    0,     0,     0,        0 }, // D.R -> D + R (km1)
				{ 0,  1,   1,  0,    0,     0,   0,  0,    0,     0,     0,        0 }, // D + RNAP -> D.RNAP (k2)
				{ 0,  0,   0,  0,    1,     0,   0,  0,    0,     0,     0,        0 }, // D.RNAP -> D + RNAP (km2)
				{ 0,  0,   0,  0,    1,     0,   0,  0,    0,     0,     0,        0 }, // D.RNAP -> TrRNAP (k3)
				{ 0,  0,   0,  0,    0,     1,   0,  0,    0,     0,     0,        0 }, // TrRNAP -> RBS + D + RNAP (k4)
				{ 0,  0,   0,  0,    0,     0,   1,  0,    0,     0,     0,        0 }, // RBS -> - (k5)
				{ 0,  0,   0,  0,    0,     0,   1,  1,    0,     0,     0,        0 }, // Rib + RBS -> Rib.RBS (k6)
				{ 0,  0,   0,  0,    0,     0,   0,  0,    1,     0,     0,        0 }, // Rib.RBS -> Rib + RBS (km6)
				{ 0,  0,   0,  0,    0,     0,   0,  0,    1,     0,     0,        0 }, // Rib.RBS -> ElRib + RBS (k7)
				{ 0,  0,   0,  0,    0,     0,   0,  0,    0,     1,     0,        0 }, // ElRib -> Rib + Protein (k8)
				{ 0,  0,   0,  0,    0,     0,   0,  0,    0,     0,     1,        0 }, // Protein -> FoldedProtein (k9)
				{ 0,  0,   0,  0,    0,     0,   0,  0,    0,     0,     1,        0 }, // Protein -> - (k10)
				{ 0,  0,   0,  0,    0,     0,   0,  0,    0,     0,     0,        1 }, // FoldedProtein -> - (k11)
		};
		int[] continuousSpecies = {};
		double N = 100;
		double deltaR = 0.5;
		double deltaS = 0.5;
		double epsilon = 0.1;
		double gamma = 0;
		double[] alpha = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
		double[] beta = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
		double t0 = 0.0;
		double t1 = 3;
		DefaultUnaryBinaryReactionNetwork net = new DefaultUnaryBinaryReactionNetwork(productionStochiometries[0].length, productionStochiometries.length);
		net.setStochiometries(productionStochiometries, consumptionStochiometries);
		net.setRateParameters(rateParameters);
		String[] speciesNames = {
				"R",
				"RNAP",
				"D",
				"D.R",
				"D.RNAP",
				"TrRNAP",
				"RBS",
				"Rib",
				"Rib.RBS",
				"ElRib",
				"Protein",
				"FoldedProtein"
		};

		this.net = net;
		this.continuousSpecies = continuousSpecies;
		this.N = N;
		this.deltaR = deltaR;
		this.deltaS = deltaS;
		this.epsilon = epsilon;
		this.gamma = gamma;
		this.alpha = alpha;
		this.beta = beta;
		this.t0 = t0;
		this.t1 = t1;
		this.x0 = x0;
		this.plotScales = plotScales;
		this.speciesNames = speciesNames;
		this.rng = new MersenneTwister();
		this.rdg = new RandomDataGenerator(this.rng);
	}

}
