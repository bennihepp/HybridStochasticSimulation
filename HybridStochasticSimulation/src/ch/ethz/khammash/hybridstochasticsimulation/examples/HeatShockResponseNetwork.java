package ch.ethz.khammash.hybridstochasticsimulation.examples;

import java.util.Arrays;

import org.apache.commons.math3.random.MersenneTwister;
import org.apache.commons.math3.random.RandomDataGenerator;

import ch.ethz.khammash.hybridstochasticsimulation.networks.DefaultUnaryBinaryReactionNetwork;


//see Kang 2012 - A multiscale approximation in a heat shock response model of E. coli
//species:
//S0: sigma32 mRNA
//S1: sigma32 protein
//S2: E sigma32
//S3: FtsH
//S4: GroEL
//S5: J
//S6: sigma32-J
//S7: Recombinant protein
//S8: J-Recombinant protein
//reactions:
//R0:  ~ -> S7             (k1) [Recombinant protein synthesis]
//R1:  S1 -> S2            (k2) [Holoenzyme association] *
//R2:  S2 -> S1            (k3) [Holoenzyme disassociation] *
//R3:  S0 -> S0 + S1       (k4) [sigma32 translation]
//R4:  S2 -> S1 + S4       (k5) [GroEL synthesis]
//R5:  S2 -> S1 + S3       (k6) [FtsH synthesis]
//R6:  S2 -> S1 + S5       (k7) [J-production] *
//R7:  S6 -> S1 + S5       (k8) [sigma32-J-disassociation]
//R8:  S1 + S5 -> S6       (k9) [sigma32-J-association]
//R9:  S5 + S7 -> S8       (k10) [Recombinant protein-J association] *
//R10: S7 -> ~             (k11) [Recombinant protein degradation]
//R11: S8 -> S5 + S7       (k12) [Recombinant protein-J disassociation] *
//R12: ~ -> S0             (k13) [sigma32 transcription]
//R13: S0 -> ~             (k14) [sigma32 mRNA decay]
//R14: S3 + S6 -> S3 + S5  (k15) [sigma32 degradation]
//R15: S4 -> ~             (k16) [GroEL degradation]
//R16: S5 -> ~             (k17) [J degradation]
//R17: S3 -> ~             (k18) [FtsH degradation]

public class HeatShockResponseNetwork extends SimulationConfiguration {

	public HeatShockResponseNetwork() {

		int[][] productionStochiometries = {
		//   S0  S1  S2  S3  S4  S5  S6  S7  S8
		    { 0,  0,  0,  0,  0,  0,  0,  1,  0,}, // R0
		    { 0,  0,  1,  0,  0,  0,  0,  0,  0,}, // R1
		    { 0,  1,  0,  0,  0,  0,  0,  0,  0,}, // R2
		    { 1,  1,  0,  0,  0,  0,  0,  0,  0,}, // R3
		    { 0,  1,  0,  0,  1,  0,  0,  0,  0,}, // R4
		    { 0,  1,  0,  1,  0,  0,  0,  0,  0,}, // R5
		    { 0,  1,  0,  0,  0,  1,  0,  0,  0,}, // R6
		    { 0,  1,  0,  0,  0,  1,  0,  0,  0,}, // R7
		    { 0,  0,  0,  0,  0,  0,  1,  0,  0,}, // R8
		    { 0,  0,  0,  0,  0,  0,  0,  0,  1,}, // R9
		    { 0,  0,  0,  0,  0,  0,  0,  0,  0,}, // R10
		    { 0,  0,  0,  0,  0,  1,  0,  1,  0,}, // R11
		    { 1,  0,  0,  0,  0,  0,  0,  0,  0,}, // R12
		    { 0,  0,  0,  0,  0,  0,  0,  0,  0,}, // R13
		    { 0,  0,  0,  1,  0,  1,  0,  0,  0,}, // R14
		    { 0,  0,  0,  0,  0,  0,  0,  0,  0,}, // R15
		    { 0,  0,  0,  0,  0,  0,  0,  0,  0,}, // R16
		    { 0,  0,  0,  0,  0,  0,  0,  0,  0,}, // R17
		};

		int[][] consumptionStochiometries = {
		//   S0  S1  S2  S3  S4  S5  S6  S7  S8
		    { 0,  0,  0,  0,  0,  0,  0,  0,  0,}, // R0
		    { 0,  1,  0,  0,  0,  0,  0,  0,  0,}, // R1
		    { 0,  0,  1,  0,  0,  0,  0,  0,  0,}, // R2
		    { 1,  0,  0,  0,  0,  0,  0,  0,  0,}, // R3
		    { 0,  0,  1,  0,  0,  0,  0,  0,  0,}, // R4
		    { 0,  0,  1,  0,  0,  0,  0,  0,  0,}, // R5
		    { 0,  0,  1,  0,  0,  0,  0,  0,  0,}, // R6
		    { 0,  0,  0,  0,  0,  0,  1,  0,  0,}, // R7
		    { 0,  1,  0,  0,  0,  1,  0,  0,  0,}, // R8
		    { 0,  0,  0,  0,  0,  1,  0,  1,  0,}, // R9
		    { 0,  0,  0,  0,  0,  0,  0,  1,  0,}, // R10
		    { 0,  0,  0,  0,  0,  0,  0,  0,  1,}, // R11
		    { 0,  0,  0,  0,  0,  0,  0,  0,  0,}, // R12
		    { 1,  0,  0,  0,  0,  0,  0,  0,  0,}, // R13
		    { 0,  0,  0,  1,  0,  0,  1,  0,  0,}, // R14
		    { 0,  0,  0,  0,  1,  0,  0,  0,  0,}, // R15
		    { 0,  0,  0,  0,  0,  1,  0,  0,  0,}, // R16
		    { 0,  0,  0,  1,  0,  0,  0,  0,  0,}, // R17
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
		    10,  // S0
		    1,   // S1
		    1,   // S2
		    93,  // S3
		    172, // S4
		    54,  // S5
		    7,   // S6
		    50,  // S7
		    0,   // S8
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
		net.setStoichiometries(productionStochiometries, consumptionStochiometries);
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
