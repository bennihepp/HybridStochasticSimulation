package ch.ethz.khammash.hybridstochasticsimulation.examples;

import java.util.Arrays;

import org.apache.commons.math3.random.MersenneTwister;
import org.apache.commons.math3.random.RandomDataGenerator;

import ch.ethz.khammash.hybridstochasticsimulation.networks.DefaultUnaryBinaryReactionNetwork;


//see Kang 2012 - A multiscale approximation in a heat shock response model of E. coli
//species:
//S1: sigma32 mRNA
//S2: sigma32 protein
//S3: E sigma32
//S4: FtsH
//S5: GroEL
//S6: J
//S7: sigma32-J
//S8: Recombinant protein
//S9: J-Recombinant protein
//reactions:
//R1:  ~ -> S8             (k1) [Recombinant protein synthesis]
//R2:  S2 -> S3            (k2) [Holoenzyme association]
//R3:  S3 -> S2            (k3) [Holoenzyme disassociation]
//R4:  S1 -> S1 + S2       (k4) [sigma32 translation]
//R5:  S3 -> S2 + S5       (k5) [GroEL synthesis]
//R6:  S3 -> S2 + S4       (k6) [FtsH synthesis]
//R7:  S3 -> S2 + S6       (k7) [J-production]
//R8:  S7 -> S2 + S6       (k8) [sigma32-J-disassociation]
//R9:  S2 + S6 -> S7       (k9) [sigma32-J-association]
//R10: S6 + S8 -> S9       (k10) [Recombinant protein-J association]
//R11: S8 -> ~             (k11) [Recombinant protein degradation]
//R12: S9 -> S6 + S8       (k12) [Recombinant protein-J disassociation]
//R13: ~ -> S1             (k13) [sigma32 transcription]
//R14: S1 -> ~             (k14) [sigma32 mRNA decay]
//R15: S4 + S7 -> S4 + S6  (k15) [sigma32 degradation]
//R16: S5 -> ~             (k16) [GroEL degradation]
//R17: S6 -> ~             (k17) [J degradation]
//R18: S4 -> ~             (k18) [FtsH degradation]

public class HeatShockResponseNetwork extends SimulationConfiguration {

	public HeatShockResponseNetwork() {

		int[][] productionStochiometries = {
			//   S1  S2  S3  S4  S5  S6  S7  S8  S9
		    { 0,  0,  0,  0,  0,  0,  0,  1,  0,}, // R1
		    { 0,  0,  1,  0,  0,  0,  0,  0,  0,}, // R2
		    { 0,  1,  0,  0,  0,  0,  0,  0,  0,}, // R3
		    { 1,  1,  0,  0,  0,  0,  0,  0,  0,}, // R4
		    { 0,  1,  0,  0,  1,  0,  0,  0,  0,}, // R5
		    { 0,  1,  0,  1,  0,  0,  0,  0,  0,}, // R6
		    { 0,  1,  0,  0,  0,  1,  0,  0,  0,}, // R7
		    { 0,  1,  0,  0,  0,  1,  0,  0,  0,}, // R8
		    { 0,  0,  0,  0,  0,  0,  1,  0,  0,}, // R9
		    { 0,  0,  0,  0,  0,  0,  0,  0,  1,}, // R10
		    { 0,  0,  0,  0,  0,  0,  0,  0,  0,}, // R11
		    { 0,  0,  0,  0,  0,  1,  0,  1,  0,}, // R12
		    { 1,  0,  0,  0,  0,  0,  0,  0,  0,}, // R13
		    { 0,  0,  0,  0,  0,  0,  0,  0,  0,}, // R14
		    { 0,  0,  0,  1,  0,  1,  0,  0,  0,}, // R15
		    { 0,  0,  0,  0,  0,  0,  0,  0,  0,}, // R16
		    { 0,  0,  0,  0,  0,  0,  0,  0,  0,}, // R17
		    { 0,  0,  0,  0,  0,  0,  0,  0,  0,}, // R18
		};

		int[][] consumptionStochiometries = {
		//   S1  S2  S3  S4  S5  S6  S7  S8  S9
		    { 0,  0,  0,  0,  0,  0,  0,  0,  0,}, // R1
		    { 0,  1,  0,  0,  0,  0,  0,  0,  0,}, // R2
		    { 0,  0,  1,  0,  0,  0,  0,  0,  0,}, // R3
		    { 1,  0,  0,  0,  0,  0,  0,  0,  0,}, // R4
		    { 0,  0,  1,  0,  0,  0,  0,  0,  0,}, // R5
		    { 0,  0,  1,  0,  0,  0,  0,  0,  0,}, // R6
		    { 0,  0,  1,  0,  0,  0,  0,  0,  0,}, // R7
		    { 0,  0,  0,  0,  0,  0,  1,  0,  0,}, // R8
		    { 0,  1,  0,  0,  0,  1,  0,  0,  0,}, // R9
		    { 0,  0,  0,  0,  0,  1,  0,  1,  0,}, // R10
		    { 0,  0,  0,  0,  0,  0,  0,  1,  0,}, // R11
		    { 0,  0,  0,  0,  0,  0,  0,  0,  1,}, // R12
		    { 0,  0,  0,  0,  0,  0,  0,  0,  0,}, // R13
		    { 1,  0,  0,  0,  0,  0,  0,  0,  0,}, // R14
		    { 0,  0,  0,  1,  0,  0,  1,  0,  0,}, // R15
		    { 0,  0,  0,  0,  1,  0,  0,  0,  0,}, // R16
		    { 0,  0,  0,  0,  0,  1,  0,  0,  0,}, // R17
		    { 0,  0,  0,  1,  0,  0,  0,  0,  0,}, // R18
		};

		double[] rateParameters = {
		    4.00e0,   // R1
		    7.00e-1,  // R2
		    1.30e-1,  // R3
		    7.00e-3,  // R4
		    6.30e-3,  // R5
		    4.88e-3,  // R6
		    4.88e-3,  // R7
		    4.40e-4,  // R8
		    3.62e-4,  // R9
		    3.62e-4,  // R10
		    9.99e-5,  // R11
		    4.40e-5,  // R12
		    1.40e-5,  // R13
		    1.40e-6,  // R14
		    1.42e-6,  // R15
		    1.80e-8,  // R16
		    6.40e-10, // R17
		    7.40e-11, // R18
		};

		double[] x0 = {
		    10,  // S1
		    1,   // S2
		    1,   // S3
		    93,  // S4
		    172, // S5
		    54,  // S6
		    7,   // S7
		    50,  // S8
		    0,   // S9
		};

		int[] importantSpecies = { 0, 1 };

		double t0 = 0.0;
		double t1 = 100.0;
		double N = 100;
//		double gamma = 0;
//		double[] alpha = new double[x0.length];
//		double[] beta = new double[rateParameters.length];
		//                R1  R2  R3  R4  R5  R6  R7  R8  R9 R10 R11 R12 R13 R14 R15 R16 R17 R18
		double[] beta = {  0,  0,  0, -1, -1, -1, -1, -2, -2, -2, -2, -2, -2, -2, -3, -2, -2, -2 };
		// Scaling 1:
//		double gamma = 0;
//		double[] alpha = { 1, 0, 0, 2, 2, 0, 0, 0, 0 };
//		double[] rho = {  0,  0,  0,  0, -1, -1, -1, -2, -2, -2, -2, -2, -2, -1, -1,  0, -2,  0 };
		// Scaling 2:
//		double gamma = 1;
//		double[] alpha = { 0, 0, 0, 2, 2, 0, 0, 1, 1 };
//		double[] rho = {  0,  0,  0, -1, -1, -1, -1, -2, -2, -1, -1, -1, -2, -2, -1,  0, -2,  0 };
		// Scaling 3:
		double gamma = 2;
		double[] alpha = { 0, 1, 1, 2, 2, 0, 0, 2, 2 };
//		double[] rho = {  0,  1,  1, -1,  0,  0,  0, -2, -1,  0,  0,  0, -2, -2, -1,  0, -2,  0 };

		DefaultUnaryBinaryReactionNetwork net = new DefaultUnaryBinaryReactionNetwork(x0.length, rateParameters.length);
		net.setStochiometries(productionStochiometries, consumptionStochiometries);
		net.setRateParameters(rateParameters);
		String[] speciesNames = {
		    "sigma32 mRNA",
		    "sigma32 protein",
		    "E sigma32",
		    "FtsH",
		    "GroEL",
		    "J",
		    "sigma32-J",
		    "Recombinant protein",
		    "J-Recombinant protein",
		};
		double[] plotScales = new double[x0.length];
		Arrays.fill(plotScales, 1.0);

		this.net = net;
		this.N = N;
		this.gamma = gamma;
		this.alpha = alpha;
		this.beta = beta;
		this.t0 = t0;
		this.t1 = t1;
		this.x0 = x0;
		this.importantSpecies = importantSpecies;
		this.plotScales = plotScales;
		this.speciesNames = speciesNames;
		this.rng = new MersenneTwister();
		this.rdg = new RandomDataGenerator(this.rng);
	}

}
