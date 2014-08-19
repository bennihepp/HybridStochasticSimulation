package ch.ethz.bhepp.hybridstochasticsimulation.gui;

import java.awt.event.ActionEvent;

public class GUIEvent {

	public static enum EventType {
		SIMULATION,
	}

	private EventType type;
	private ActionEvent event;

	public GUIEvent(EventType type) {
		this(type, null);
	}

	public GUIEvent(EventType type, ActionEvent event) {
		this.type = type;
		this.event = event;
	}

	public EventType getDescription() {
		return type;
	}

	public ActionEvent getEvent() {
		return event;
	}

}