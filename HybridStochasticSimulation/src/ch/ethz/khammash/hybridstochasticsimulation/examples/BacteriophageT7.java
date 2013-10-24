package ch.ethz.khammash.hybridstochasticsimulation.examples;

import java.util.Arrays;

import javax.inject.Inject;
import javax.inject.Named;

import org.apache.commons.math3.random.MersenneTwister;
import org.apache.commons.math3.random.RandomDataGenerator;

import ch.ethz.khammash.hybridstochasticsimulation.networks.DefaultUnaryBinaryReactionNetwork;


// Intracellular growth of bacteriophage T7.
// Original publication: Srivastava et al., "Stochastic vs. deterministic modeling of intracellular viral kinetics.", J. theor Biol, 218:309-321 (2002).
//                  and  Alfonsi et al., "Adaptive simulation of hybrid stochastic and deterministic models for biochemical systems.", ESAIM: Proceedings, Vol. 14 (2005).
// species:
//  S0: tem
//  S1: gen
//  S2: struc
//  S3: virus
// reactions:
//  R0:  gen         -> tem            [0.025]
//  R1:  tem         -> -              [0.25]
//  R2:  tem         -> tem + gen      [1.0]
//  R3:  gen + struc -> virus          [7.5e-6]
//  R4:  tem         -> tem + struc    [1000.0]
//  R5:  struc       -> -              [1.99]

public class BacteriophageT7 extends SimulationConfiguration {

	@Inject @Named("BacteriophageT7")
	public BacteriophageT7() {

		int[][] productionStochiometries = {
		//   R0  R1  R2  R3  R4  R5
		    { 1,  0,  1,  0,  1,  0 }, // S0: tem
		    { 0,  0,  1,  0,  0,  0 }, // S1: gen
		    { 0,  0,  0,  0,  1,  0 }, // S2: struc
		    { 0,  0,  0,  1,  0,  0 }, // S3: virus
		};

		int[][] consumptionStochiometries = {
		//   R0  R1  R2  R3  R4  R5
		    { 0,  1,  1,  0,  1,  0 }, // S0: tem
		    { 1,  0,  0,  1,  0,  0 }, // S1: gen
		    { 0,  0,  0,  1,  0,  1 }, // S2: struc
		    { 0,  0,  0,  0,  0,  0 }, // S3: virus
		};

		double[] x0 = { 1, 0, 0, 0 };

		double[] rateParameters = new double[productionStochiometries[0].length];
		/* Reaction propensities */
		rateParameters[0] = 0.025;
		rateParameters[1] = 0.25;
		rateParameters[2] = 1.0;
		rateParameters[3] = 7.5e-6;
		rateParameters[4] = 1000.0;
		rateParameters[5] = 1.99;

		double t0 = 0.0;
		double t1 = 200.0;
		double N = 100;
		double gamma = 0;
		double[] alpha = new double[x0.length];
		double[] beta = new double[rateParameters.length];
		int[] importantSpecies = {0};

		productionStochiometries = transpose(productionStochiometries);
		consumptionStochiometries = transpose(consumptionStochiometries);
		DefaultUnaryBinaryReactionNetwork net = new DefaultUnaryBinaryReactionNetwork(x0.length, rateParameters.length);
		net.setStoichiometries(productionStochiometries, consumptionStochiometries);
		net.setRateParameters(rateParameters);
		String[] speciesNames = {
		    "tem",
		    "gen",
		    "struc",
		    "virus",
		};
		double[] plotScales = new double[x0.length];
		Arrays.fill(plotScales, 1.0);

		this.net = net;
		this.N = N;
		this.gamma = gamma;
		this.alpha = alpha;
		this.beta = beta;
		this.importantSpecies = importantSpecies;
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
