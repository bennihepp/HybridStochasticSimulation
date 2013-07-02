package ch.ethz.khammash.hybridstochasticsimulation.examples;

import java.util.Arrays;

import org.apache.commons.math3.random.MersenneTwister;
import org.apache.commons.math3.random.RandomDataGenerator;

import ch.ethz.khammash.hybridstochasticsimulation.networks.DefaultUnaryBinaryReactionNetwork;


// Chemical Oscillator example.
// Original publication: Vilar, Kueh, Barkai and Liebler, "Mechanisms of noise-resistance in genetic oscillators." PNAS, Vol 99 (2002).
// species:
//  S0: A
//  S1: A_R
//  S2: Pa
//  S3: R
//  S4: Pa_A
//  S5: Pr
//  S6: Pr_A
//  S7: mRNA_a
//  S8: mRNA_r
// reactions:
//  R0:  Pa     -> Pa     + mRNA_a     [50.0]
//  R1:  Pa_A   -> Pa_A   + mRNA_a     [500.0]
//  R2:  Pr     -> Pr     + mRNA_r     [0.01]
//  R3:  Pr_A   -> Pr_A   + mRNA_r     [50.0]
//  R4:  mRNA_a -> mRNA_a + A          [50.0]
//  R5:  mRNA_r -> mRNA_r + R          [5.0]
//  R6:  A + R  -> A_R                 [20.0]
//  R7:  A + Pa -> Pa_A                [1.0]
//  R8:  Pa_A   -> A      + Pa         [50.0]
//  R9:  A + Pr -> Pr_A                [1.0]
//  R10: Pr_A   -> A      + Pr         [100.0]
//  R11: A      -> ~                   [1.0]
//  R12: R      -> ~                   [0.2]
//  R13: mRNA_a -> ~                   [10.0]
//  R14: mRNA_r -> ~                   [0.5]
//  R15: A_R    -> R                   [1.0]

public class VilarOscillator extends ExampleConfiguration {

	public VilarOscillator() {
		this(false);
	}

	public VilarOscillator(boolean modifiedParameters) {

//		int i_A = 0;
//		int i_A_R = 1;
		int i_Pa = 2;
//		int i_R = 3;
//		int i_Pa_A = 4;
		int i_Pr = 5;
//		int i_Pr_A = 6;
//		int i_mRNA_a = 7;
//		int i_mRNA_r = 8;

		int[][] productionStochiometries = {
		//   R0  R1  R2  R3  R4  R5  R6  R7  R8  R9 R10 R11 R12 R13 R14 R15
		    { 0,  0,  0,  0,  1,  0,  0,  0,  1,  0,  1,  0,  0,  0,  0,  0 }, // S0: A
		    { 0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0 }, // S1: A_R
		    { 1,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0 }, // S2: Pa
		    { 0,  0,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1 }, // S3: R
		    { 0,  1,  0,  0,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0 }, // S4: Pa_A
		    { 0,  0,  1,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  0,  0 }, // S5: Pr
		    { 0,  0,  0,  1,  0,  0,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0 }, // S6: Pr_A
		    { 1,  1,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0 }, // S7: mRNA_a
		    { 0,  0,  1,  1,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0 }, // S8: mRNA_r
		};

		int[][] consumptionStochiometries = {
		//   R0  R1  R2  R3  R4  R5  R6  R7  R8  R9 R10 R11 R12 R13 R14 R15
		    { 0,  0,  0,  0,  0,  0,  1,  1,  0,  1,  0,  1,  0,  0,  0,  0 }, // S0: A
		    { 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1 }, // S1: A_R
		    { 1,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0 }, // S2: Pa
		    { 0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  0,  0,  1,  0,  0,  0 }, // S3: R
		    { 0,  1,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0 }, // S4: Pa_A
		    { 0,  0,  1,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0 }, // S5: Pr
		    { 0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  0,  0 }, // S6: Pr_A
		    { 0,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0 }, // S7: mRNA_a
		    { 0,  0,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0 }, // S8: mRNA_r
		};

		double[] x0 = new double[productionStochiometries.length];
		x0[i_Pa] = 1;
		x0[i_Pr] = 1;

		double[] rateParameters = new double[productionStochiometries[0].length];
		// Original parameters
		double alphaA = 50.0;
		double alphaR = 0.01;
		double betaA = 50.0;
		double betaR = 5.0;
		double gammaC = 20.0;
		double gammaA = 1.0;
		double thetaA = 50.0;
		double gammaR = 1.0;
		double thetaR = 100.0;
		double deltaA = 1.0;
		double deltaR = 0.20;
		double deltaMA = 10.0;
		double deltaMR = 0.5;
		double delta_C = 1.0;
		double alpha_a = 10.0;
		double alpha_r = 5000.0;
		// Modified parameters
		if (modifiedParameters) {
			alpha_r = 50;
			gammaC = 0.2;
			alphaA *= 10.0;
			alphaR *= 10.0;
//			deltaMA *= 10.0;
//			deltaMR *= 10.0;
//	
//	//		alphaA = 5.0;
//	//		alpha_a = 100.0;
//	//		deltaMA = 1.0;
//			betaA = 100.0;
//			gammaR = 0.1;
//			delta_C = 0.1;
//	//		gammaR = 0.1;
//	//		thetaR = 10.0;
		}
		// Reaction propensities
		rateParameters[0] = alphaA;            // 50.0
		rateParameters[1] = alpha_a * alphaA;  // 500.0
		rateParameters[2] = alphaR;            // 0.01
		rateParameters[3] = alpha_r * alphaR;  // 50.0
		rateParameters[4] = betaA;             // 50.0
		rateParameters[5] = betaR;             // 5.0
		rateParameters[6] = gammaC;            // 20.0
		rateParameters[7] = gammaA;            // 1.0
		rateParameters[8] = thetaA;            // 50.0
		rateParameters[9] = gammaR;            // 1.0
		rateParameters[10] = thetaR;           // 100.0
		rateParameters[11] = deltaA;           // 1.0
		rateParameters[12] = deltaR;           // 0.2
		rateParameters[13] = deltaMA;          // 10.0
		rateParameters[14] = deltaMR;          // 0.5
		rateParameters[15] = delta_C;          // 1.0

		double t0 = 0.0;
		double t1 = 100.0;
		double N = 100;
		double gamma = 0;
		double[] alpha = new double[x0.length];
		double[] beta = new double[rateParameters.length];
//		int[] importantSpecies = { 0 };

		productionStochiometries = transpose(productionStochiometries);
		consumptionStochiometries = transpose(consumptionStochiometries);
		DefaultUnaryBinaryReactionNetwork net = new DefaultUnaryBinaryReactionNetwork(x0.length, rateParameters.length);
		net.setStochiometries(productionStochiometries, consumptionStochiometries);
		net.setRateParameters(rateParameters);
		String[] speciesNames = {
		    "A",
		    "A_R",
		    "Pa",
		    "R",
		    "Pa_A",
		    "Pr",
		    "Pr_A",
		    "mRNA_a",
		    "mRNA_r",
		};
		double[] plotScales = new double[x0.length];
		Arrays.fill(plotScales, 1.0);

		this.net = net;
		this.N = N;
		this.gamma = gamma;
		this.alpha = alpha;
		this.beta = beta;
//		this.importantSpecies = importantSpecies;
		this.t0 = t0;
		this.t1 = t1;
		this.x0 = x0;
		this.plotScales = plotScales;
		this.speciesNames = speciesNames;
		this.rng = new MersenneTwister();
		this.rdg = new RandomDataGenerator(this.rng);
	}

	private int[][] transpose(int[][] matrix) {
		int[][] result = new int[matrix[0].length][matrix.length];
		for (int i=0; i < matrix.length; i++)
			for (int j=0; j < matrix[i].length; j++)
				result[j][i] = matrix[i][j];
		return result;
	}

}
