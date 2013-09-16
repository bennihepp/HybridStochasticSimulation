import ch.ethz.khammash.hybridstochasticsimulation.models.ReactionNetworkModel
import ch.ethz.khammash.hybridstochasticsimulation.batch.SimulationJob
import ch.ethz.khammash.hybridstochasticsimulation.trajectories.FiniteTrajectory

class ScalaSimulationJob[T >: ReactionNetworkModel] extends SimulationJob {

	def getRuns() = 5

	def runJob() {}

	def runSingleSimulation() = null

	def beginOutput() {}

	def endOutput() {}

	def addSimulationResult(tr: FiniteTrajectory) {}

}
