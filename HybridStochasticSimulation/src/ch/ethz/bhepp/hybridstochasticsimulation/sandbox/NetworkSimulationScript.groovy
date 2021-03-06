package ch.ethz.khammash.hybridstochasticsimulation.sandbox

import ch.ethz.khammash.hybridstochasticsimulation.GUI;
import ch.ethz.khammash.hybridstochasticsimulation.Utilities;
import ch.ethz.khammash.hybridstochasticsimulation.gui.PlotWindow

class NetworkSimulationScript {

	static main(args) {
		def plots = regulatedTranscriptionNetwork()
		def rows = (int)(plots.size() / 2)
		def cols = (int)(plots.size() / rows)
		def window = new PlotWindow("PDMP simulations", rows, cols)
		for (plot in plots)
			window.add(plot)
		window.setVisible(true)
	}

	static regulatedTranscriptionNetwork() {
		//def nss = NetworkSimulation.loadSimpleCrystallizationNetwork()
		def nss = GUI.loadRegulatedTranscriptionNetwork()
		
//		nss.rng = new MersenneTwister(100)
//		nss.rdg = new RandomDataGenerator(nss.rng)

		def plots = [] 

		def runs = 100
		def numberOfTimePoints = 1001

		def tVector = Utilities.computeTimeVector(numberOfTimePoints,
			nss.t0, nss.t1)

		def std = GUI.simulateStochasticDistribution(
			runs, nss, tVector)
		def plot = GUI.plotTrajectoryDistribution(nss, std,)
		plot.setTitle("Stochastic")
		plots.add(plot)

		def st = GUI.simulateStochastic(nss)
		plot = GUI.plotTrajectory(nss, st)
		plot.setTitle("Stochastic single trajectory")
		plots.add(plot)

		std = GUI.simulateMSPDMPDistribution(runs, nss, tVector)
		plot = GUI.plotTrajectoryDistribution(nss, std)
		plot.setTitle("PDMP")
		plots.add(plot)

		st = GUI.simulateMSPDMP(nss, tVector)
		plot = GUI.plotTrajectory(nss, st)
		plot.setTitle("PDMP single trajectory")
		plots.add(plot)

		nss.gamma = 1
		nss.t1 = 1000.0
		
		tVector = Utilities.computeTimeVector(numberOfTimePoints,
			nss.t0, nss.t1)

		std = GUI.simulateStochasticDistribution(
			runs, nss, tVector)
		plot = GUI.plotTrajectoryDistribution(nss, std)
		plot.setTitle("Stochatic, gamma=1")
		plots.add(plot)

		std = GUI.simulateMSPDMPDistribution(runs, nss, tVector)
		plot = GUI.plotTrajectoryDistribution(nss, std)
		plot.setTitle("PDMP, gamma=1")
		plots.add(plot)

		return plots
	}

}
