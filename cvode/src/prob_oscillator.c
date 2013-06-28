/*
 *  prob_oscillator.c
 *  
 *	Includes the stoichiometries and propensities for the problem.
 * 	Chemical Oscillator example
 *      9 chemical species
 *	16 chemical reactions
 *	16 rate parameters
 *  Original reference:
 *      Vilar, Kueh, Barkai and Liebler, "Mechanisms of noise-resistance  
 *			in genetic oscillators." PNAS, Vol 99 (2002).
 * Used as an example in the following references:
 * M. Rathinam, P. Sheppard, and M. Khammash, J. Chem. Phys. 132, 034103 (2010).
 * P. Sheppard, M. Rathinam, and M. Khammash, J. Chem. Phys. 136, 034115 (2012).

 SPSens : Stochastic parameter sensitivity analysis for chemical networks
 Copyright (C) 2012  Patrick W. Sheppard
 
 */
 
const int N_SPECIES = 9;
const int N_RXNS = 16;
const int N_PARAMS = 16;

void initialize_state (double* X)
{
	const int A = 0;
	const int A_R = 1;
	const int Pa = 2;
	const int R = 3;
	const int Pa_A = 4;
	const int Pr = 5;
	const int Pr_A = 6;
	const int mRNA_a = 7;
	const int mRNA_r = 8;
	
	X[A] = 0;
	X[A_R] = 0;
	X[Pa] = 1;
	X[R] = 0;
	X[Pa_A] = 0;
	X[Pr] = 1;
	X[Pr_A] = 0;
	X[mRNA_a] = 0; 
	X[mRNA_r] = 0;
}

void initialize_parameters(double* c)
{
	/* Parameters */
	double alphaA;
	double alphaR;
	double betaA;
	double betaR;
	double gammaC;
	double betaC;
	double gammaA;
	double thetaA;
	double gammaR;
	double thetaR;
	double deltaA; 
	double deltaR;
	double deltaMA;
	double deltaMR;
	double delta_C;
	double alpha_a;
	double alpha_r;
	
	c[0] = alphaA = 50.0;
	c[1] = alphaR = 0.01;
	c[2] = betaA = 50.0;
	c[3] = betaR = 5.0;
	c[4] = gammaC = 20.0;
	c[5] = gammaA = 1.0;
	c[6] = thetaA = 50.0;
	c[7] = gammaR = 1.0;
	c[8] = thetaR = 100.0;
	c[9] = deltaA = 1.0;
	c[10] = deltaR = 0.20;
	c[11] = deltaMA = 10.0;
	c[12] = deltaMR = 0.5;
	c[13] = delta_C = 1.0;
	c[14] = alpha_a = 10.0;
	c[15] = alpha_r = 5000.0;
}


void state_change(int* delta_x, const int rxn)
{
	int idx_s;	
	for (idx_s = 0; idx_s<N_SPECIES; idx_s++){
		delta_x[idx_s]=0;
	}
	const int A = 0;
	const int A_R = 1;
	const int Pa = 2;
	const int R = 3;
	const int Pa_A = 4;
	const int Pr = 5;
	const int Pr_A = 6;
	const int mRNA_a = 7;
	const int mRNA_r = 8;
	if (rxn==0){
		delta_x[mRNA_a] = 1;
	} else if (rxn==1)	{
		delta_x[mRNA_a] = 1;
	} else if (rxn==2)	{
		delta_x[mRNA_r] = 1;
	} else if (rxn==3)	{
		delta_x[mRNA_r] = 1;
	} else if (rxn==4)	{
		delta_x[A] = 1;
	} else if (rxn==5)	{
		delta_x[R] = 1;
	} else if (rxn==6)	{
		delta_x[A] = -1;
		delta_x[R] = -1;
		delta_x[A_R] = 1;
	} else if (rxn==7)	{
		delta_x[A] = -1;
		delta_x[Pa] = -1;
		delta_x[Pa_A] = 1;
	} else if (rxn==8)	{
		delta_x[A] = 1;
		delta_x[Pa] = 1;
		delta_x[Pa_A] = -1;
	} else if (rxn==9)	{
		delta_x[A] = -1;
		delta_x[Pr] = -1;
		delta_x[Pr_A] = 1;
	} else if (rxn==10)	{
		delta_x[A] = 1;
		delta_x[Pr] = 1;
		delta_x[Pr_A] = -1;
	} else if (rxn==11)	{
		delta_x[A] = -1;
	} else if (rxn==12)	{
		delta_x[R] = -1;
	} else if (rxn==13)	{
		delta_x[mRNA_a] = -1;
	} else if (rxn==14)	{
		delta_x[mRNA_r] = -1;
	} else if (rxn==15)	{
		delta_x[A_R] = -1;
		delta_x[R] = 1;
	}
}
void propensity(const double* X, const double* params, double* prop)	
{
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
	/* Species */
	double A = (double) X[0];
	double A_R = (double) X[1];
	double Pa = (double) X[2];
	double R = (double) X[3];
	double Pa_A = (double) X[4];
	double Pr = (double) X[5];
	double Pr_A = (double) X[6];
	double mRNA_a = (double) X[7];
	double mRNA_r = (double) X[8];
	/* Reaction propensities */
	prop[0] = alphaA*Pa;
	prop[1] = alpha_a*alphaA*Pa_A;		
	prop[2] = alphaR*Pr;
	prop[3] = alpha_r*alphaR*Pr_A;
	prop[4] = betaA*mRNA_a;
	prop[5] = betaR*mRNA_r;
	prop[6] = gammaC * A * R;
	prop[7] = gammaA * Pa * A;
	prop[8] = thetaA*Pa_A;
	prop[9] = gammaR * Pr * A;
	prop[10] = thetaR*Pr_A;
	prop[11] = deltaA*A;
	prop[12] = deltaR*R;
	prop[13] = deltaMA*mRNA_a;
	prop[14] = deltaMR*mRNA_r;
	prop[15] = delta_C*A_R;
    
}



const int N_OUTPUTS=3; 
int output_function(const double t, const int* x, const double *c, double* output, int count)	
{
	const int A = 0;
	const int A_R = 1;
	const int Pa = 2;
	const int R = 3;
	const int Pa_A = 4;
	const int Pr = 5;
	const int Pr_A = 6;
	const int mRNA_a = 7;
	const int mRNA_r = 8;
	/* In this example, return all activity regarding the repressor */
	output[0] = (double) x[R];
	output[1] = (double) x[A_R];
	output[2] = (double) x[mRNA_r];
	return 0;
}
