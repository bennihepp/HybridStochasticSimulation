package ch.ethz.khammash.hybridstochasticsimulation.controllers;

import ch.ethz.khammash.hybridstochasticsimulation.models.PDMPModel;

public interface PDMPModelFactory {
	public PDMPModel createModel();
}