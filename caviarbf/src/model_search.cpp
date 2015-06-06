#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>

#include "model_selection.h"
#include "bf_io.h"
#include "InputArgument.h"

using namespace std;

int main(int argc, char** argv) {
	InputArgumentModelSearch inputArgument;
	inputArgument.processArguments(argc, argv);
	cout << inputArgument.priorFile;

	double rho = 1;
	
	string marginalFile = inputArgument.outputFile + ".marginal";
	ofstream marginalOutput(marginalFile.c_str());
	string statisticsFile = inputArgument.outputFile + ".statistics";
	ofstream statisticsOutput(statisticsFile.c_str());

//	string m1stepwiseFile = string(outputFile) + ".m.stepwise";
//	ofstream m1stepwiseOutput(m1stepwiseFile.c_str());

	string logFile = inputArgument.inputFile + ".log";

	int max_n = 1000;
	int max_k = 10;
	CombinationTable combinationTable(max_n,  max_k);
	try {
		vector<double> allSNPPriors;
		if (inputArgument.priorInAFile) {
			readPriorFile(inputArgument.priorFile.c_str(), \
			              inputArgument.nSNPs, allSNPPriors);
		}
		BayesFactorData bfData(inputArgument.inputFile.c_str(), \
		                       inputArgument.nSNPs);
		double* modelPrior = new double[bfData.nModel];
		double* modelScore = new double[bfData.nModel];
		double nullModelPrior;
		if (inputArgument.priorInAFile) {
			calculatePriors(bfData, allSNPPriors, modelPrior, nullModelPrior);
		} else {
			calculatePriors(bfData, inputArgument.prior, \
			                modelPrior, nullModelPrior);
		}
		double totalScore;
		calculateModelScoreAndSum(bfData, modelPrior, modelScore, \
									 nullModelPrior, totalScore);
		outputStatistics(totalScore, nullModelPrior, statisticsOutput);
		vector<int> selected; // = SNPList.models.back().SNPSet;

		StepwiseSelection marginalProbabilitySelection;
		marginalProbability(bfData, modelScore, totalScore, \
		                    marginalProbabilitySelection);
		marginalOutput << marginalProbabilitySelection;
		
		if (inputArgument.outputStepwiseResult) {
			string stepwiseFile = inputArgument.outputFile + ".stepwise";
			ofstream stepwiseOutput(stepwiseFile.c_str());
			StepwiseSelection stepwiseResult;
			greedySearch(bfData, modelScore, totalScore, \
						 selected, rho, stepwiseResult);
			stepwiseOutput << stepwiseResult;
		}
		
		if (inputArgument.outputExhaustiveResult || \
			inputArgument.outputExhaustiveStepwiseResult) {
			string exhaustiveFile = inputArgument.outputFile + ".exhaustive";
			ofstream exhaustiveOutput(exhaustiveFile.c_str());
			ExhaustiveSelection SNPList;
			exhaustiveSearch(bfData,  modelScore, totalScore, \
							 bfData.maxInModel, SNPList);
			exhaustiveOutput << SNPList;
			
			if (inputArgument.outputExhaustiveStepwiseResult) {
				string exhaustivestepwiseFile = inputArgument.outputFile + \
									".exhaustivestepwise";
				ofstream exhaustivestepwiseOutput( \
										exhaustivestepwiseFile.c_str());
				selected = SNPList.models.back().SNPSet;
				StepwiseSelection exhaustivestepwiseResult;
				greedySearch(bfData, modelScore, totalScore, \
							 selected, rho, exhaustivestepwiseResult);
				ExhaustiveSelection SNPList2;
				combineExhaustiveStepwise(SNPList, exhaustivestepwiseResult, \
											 SNPList2);
				exhaustivestepwiseOutput << SNPList2;
			}
		}
		
//		selected = marginalProbabilitySelection.selectedOfEachStep;
//		selected.resize(1); // keep the first
//		StepwiseSelection m1stepwiseResult;
//		greedySearch(bfData, modelScore, totalScore, \
//					 selected, rho, m1stepwiseResult);
//		m1stepwiseOutput << m1stepwiseResult;

		delete [] modelPrior;
		delete [] modelScore;
	} catch (std::exception& e) {
		std::cerr << e.what() << endl;
		return 1;
	}
	return 0;
}