package ch.ethz.khammash.hybridstochasticsimulation

import ch.ethz.khammash.hybridstochasticsimulation.examples.NetworkSimulation;
import ch.ethz.khammash.hybridstochasticsimulation.plotting.GridWindow;

class NetworkSimulationScript {

	static main(args) {
		def plots = regulatedTranscriptionNetwork()
		def rows = (int)(plots.size() / 2)
		def cols = (int)(plots.size() / rows)
		def window = new GridWindow("PDMP simulations", rows, cols)
		for (plot in plots)
			window.add(plot)
		window.setVisible(true)
	}

	static regulatedTranscriptionNetwork() {
		//def nss = NetworkSimulation.loadSimpleCrystallizationNetwork()
		def nss = NetworkSimulation.loadRegulatedTranscriptionNetwork()
		
//		nss.rng = new MersenneTwister(100)
//		nss.rdg = new RandomDataGenerator(nss.rng)

		def plots = [] 

		def runs = 100
		def numberOfTimePoints = 1001

		def tVector = Utilities.computeTimeVector(numberOfTimePoints,
			nss.t0, nss.t1)

		def std = NetworkSimulation.simulateStochasticDistribution(
			runs, nss, tVector)
		def plot = NetworkSimulation.plotTrajectoryDistribution(nss, std,)
		plot.setTitle("Stochastic")
		plots.add(plot)

		def st = NetworkSimulation.simulateStochastic(nss)
		plot = NetworkSimulation.plotTrajectory(nss, st)
		plot.setTitle("Stochastic single trajectory")
		plots.add(plot)

		std = NetworkSimulation.simulateMSPDMPDistribution(runs, nss, tVector)
		plot = NetworkSimulation.plotTrajectoryDistribution(nss, std)
		plot.setTitle("PDMP")
		plots.add(plot)

		st = NetworkSimulation.simulateMSPDMP(nss, tVector)
		plot = NetworkSimulation.plotTrajectory(nss, st)
		plot.setTitle("PDMP single trajectory")
		plots.add(plot)

		nss.gamma = 1
		nss.t1 = 1000.0
		
		tVector = Utilities.computeTimeVector(numberOfTimePoints,
			nss.t0, nss.t1)

		std = NetworkSimulation.simulateStochasticDistribution(
			runs, nss, tVector)
		plot = NetworkSimulation.plotTrajectoryDistribution(nss, std)
		plot.setTitle("Stochatic, gamma=1")
		plots.add(plot)

		std = NetworkSimulation.simulateMSPDMPDistribution(runs, nss, tVector)
		plot = NetworkSimulation.plotTrajectoryDistribution(nss, std)
		plot.setTitle("PDMP, gamma=1")
		plots.add(plot)

		return plots
	}

}
