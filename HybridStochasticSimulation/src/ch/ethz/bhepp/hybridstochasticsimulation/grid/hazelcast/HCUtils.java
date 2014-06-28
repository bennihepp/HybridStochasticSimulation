package ch.ethz.bhepp.hybridstochasticsimulation.grid.hazelcast;

import java.io.FileNotFoundException;

import ch.ethz.bhepp.hybridstochasticsimulation.grid.hazelcast.HCController.Control;
import ch.ethz.bhepp.hybridstochasticsimulation.trajectories.FiniteTrajectory;

import com.hazelcast.config.FileSystemXmlConfig;
import com.hazelcast.core.Hazelcast;
import com.hazelcast.core.HazelcastInstance;
import com.hazelcast.core.IQueue;
import com.hazelcast.core.ITopic;

public class HCUtils {

	private static final String DEFAULT_HAZELCAST_CONFIG_FILE = "hazelcast.xml";

	private static final String TASK_QUEUE_NAME = "tasks";

	private static final String RESULT_QUEUE_NAME = "results";

	private static final String CONTROL_TOPIC_NAME = "controls";

	public static HCUtils createInstance() throws FileNotFoundException {
		return createInstance(null);
	}

	public static HCUtils createInstance(String hazelcastConfigFilename) throws FileNotFoundException {
		if (hazelcastConfigFilename == null)
			hazelcastConfigFilename = DEFAULT_HAZELCAST_CONFIG_FILE;
		return new HCUtils(hazelcastConfigFilename);
	}

	private FileSystemXmlConfig hazelcastConfig;
	private volatile HazelcastInstance hazelcast;

	private HCUtils(String hazelcastConfigFilename) throws FileNotFoundException {
		this.hazelcastConfig = new FileSystemXmlConfig(hazelcastConfigFilename);
    }

	public HazelcastInstance getHazelcastInstance() {
		if (hazelcast == null) {
			synchronized (this) {
				if (hazelcast == null) {
					hazelcast = Hazelcast.newHazelcastInstance(hazelcastConfig);
				}
			}
		}
		return hazelcast;
	}

	public IQueue<Integer> getTaskQueue() {
		IQueue<Integer> taskQueue = getHazelcastInstance().getQueue(TASK_QUEUE_NAME);
		return taskQueue;
//		return getHazelcastInstance().getQueue(TASK_QUEUE_NAME);
	}

	public IQueue<FiniteTrajectory> getResultQueue() {
		return getHazelcastInstance().getQueue(RESULT_QUEUE_NAME);
	}

	public ITopic<Control> getControlTopic() {
		return getHazelcastInstance().getTopic(CONTROL_TOPIC_NAME);
	}

	public void shutdown() {
		getHazelcastInstance().getLifecycleService().shutdown();
	}

}
