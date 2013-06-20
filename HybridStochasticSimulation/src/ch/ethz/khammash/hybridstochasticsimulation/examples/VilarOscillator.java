package ch.ethz.khammash.hybridstochasticsimulation.examples;

import java.util.Arrays;

import org.apache.commons.math3.random.MersenneTwister;
import org.apache.commons.math3.random.RandomDataGenerator;

import ch.ethz.khammash.hybridstochasticsimulation.networks.DefaultUnaryBinaryReactionNetwork;


// Chemical Oscillator example.
// Original publication: Vilar, Kueh, Barkai and Liebler, "Mechanisms of noise-resistance in genetic oscillators." PNAS, Vol 99 (2002).
// species:
//  S1: A
//  S2: A_R
//  S3: Pa
//  S4: R
//  S5: Pa_A
//  S6: Pr
//  S7: Pr_A
//  S8: mRNA_a
//  S9: mRNA_r
// reactions:
//  R1:  Pa     -> Pa     + mRNA_a
//  R2:  Pa_A   -> Pa_A   + mRNA_a
//  R3:  Pr     -> Pr     + mRNA_r
//  R4:  Pr_A   -> Pr_A   + mRNA_r
//  R5:  mRNA_a -> mRNA_a + A        *
//  R6:  mRNA_r -> mRNA_r + R
//  R7:  A + R  -> A_R               *
//  R8:  A + Pa -> Pa_A
//  R9:  Pa_A   -> A      + Pa
//  R10: A + Pr -> Pr_A
//  R11: Pr_A   -> A      + Pr
//  R12: A      -> ~
//  R13: R      -> ~
//  R14: mRNA_a -> ~
//  R15: mRNA_r -> ~
//  R16: A_R    -> R                 *

public class VilarOscillator extends ExampleConfiguration {

	public VilarOscillator() {

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
		//   R1  R2  R3  R4  R5  R6  R7  R8  R9 R10 R11 R12 R13 R14 R15 R16
		    { 0,  0,  0,  0,  1,  0,  0,  0,  1,  0,  1,  0,  0,  0,  0,  0 }, // S1: A
		    { 0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0 }, // S2: A_R
		    { 1,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0 }, // S3: Pa
		    { 0,  0,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1 }, // S4: R
		    { 0,  1,  0,  0,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0 }, // S5: Pa_A
		    { 0,  0,  1,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  0,  0 }, // S6: Pr
		    { 0,  0,  0,  1,  0,  0,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0 }, // S7: Pr_A
		    { 1,  1,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0 }, // S8: mRNA_a
		    { 0,  0,  1,  1,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0 }, // S9: mRNA_r
		};

		int[][] consumptionStochiometries = {
		//   R1  R2  R3  R4  R5  R6  R7  R8  R9 R10 R11 R12 R13 R14 R15 R16
		    { 0,  0,  0,  0,  0,  0,  1,  1,  0,  1,  0,  1,  0,  0,  0,  0 }, // S1: A
		    { 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1 }, // S2: A_R
		    { 1,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0 }, // S3: Pa
		    { 0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  0,  0,  1,  0,  0,  0 }, // S4: R
		    { 0,  1,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0 }, // S5: Pa_A
		    { 0,  0,  1,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0 }, // S6: Pr
		    { 0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  0,  0 }, // S7: Pr_A
		    { 0,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0 }, // S8: mRNA_a
		    { 0,  0,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0 }, // S9: mRNA_r
		};

		double[] x0 = new double[productionStochiometries.length];
		x0[i_Pa] = 1;
		x0[i_Pr] = 1;

		double[] rateParameters = new double[productionStochiometries[0].length];
		double[] params = new double[16];
		params[0] = 50.0;
		params[1] = 0.01;
		params[2] = 50.0;
		params[3] = 5.0;
		params[4] = 20.0;
		params[5] = 1.0;
		params[6] = 50.0;
		params[7] = 1.0;
		params[8] = 100.0;
		params[9] = 1.0;
		params[10] = 0.20;
		params[11] = 10.0;
		params[12] = 0.5;
		params[13] = 1.0;
		params[14] = 10.0;
		params[15] = 5000.0;
		/* Parameters */
		double alphaA = params[0];
		double alphaR = params[1];
		double betaA = params[2];
		double betaR = params[3];
		double gammaC = params[4];
		double gammaA = params[5];
		double thetaA = params[6];
		double gammaR = params[7];
		double thetaR = params[8];
		double deltaA = params[9];
		double deltaR = params[10];
		double deltaMA = params[11];
		double deltaMR = params[12];
		double delta_C = params[13];
		double alpha_a = params[14];
		double alpha_r = params[15];
		/* Reaction propensities */
		rateParameters[0] = alphaA;
		rateParameters[1] = alpha_a*alphaA;		
		rateParameters[2] = alphaR;
		rateParameters[3] = alpha_r*alphaR;
		rateParameters[4] = betaA;
		rateParameters[5] = betaR;
		rateParameters[6] = gammaC;
		rateParameters[7] = gammaA;
		rateParameters[8] = thetaA;
		rateParameters[9] = gammaR;
		rateParameters[10] = thetaR;
		rateParameters[11] = deltaA;
		rateParameters[12] = deltaR;
		rateParameters[13] = deltaMA;
		rateParameters[14] = deltaMR;
		rateParameters[15] = delta_C;

		double t0 = 0.0;
		double t1 = 100.0;
		double N = 100;
		double gamma = 0;
		double[] alpha = new double[x0.length];
		double[] beta = new double[rateParameters.length];

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
