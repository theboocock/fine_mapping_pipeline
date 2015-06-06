#include <unistd.h>

#include "InputArgument.h"

using namespace std;

const string version = "0.1";

InputArgumentModelSearch::InputArgumentModelSearch() {
}

InputArgumentModelSearch::~InputArgumentModelSearch()
{
}

void InputArgumentModelSearch::processArguments(int argc, char** argv) {
	try {  
		TCLAP::CmdLine cmd("Search models and output probabilities "
		                   "based on Bayes factors", ' ', version);
		TCLAP::ValueArg<std::string> outputFileArg("o","output", \
		                        "output file prefix", \
								    true, "", "string");
		TCLAP::ValueArg<std::string> inputFileArg("i","input", \
								"input file storing Bayes factors", \
								true, "", "string");
		TCLAP::ValueArg<int> nSNPsArg("m","snp-number", \
		                        "the total number of variants in the data", \
								    true, 0, "integer");
		TCLAP::ValueArg<double> priorArg("p","prior", \
		"the prior probability of each SNP being causal, in the range [0, 1)\n" 
		 "if it is 0, the prior probability will be set to 1 / m, \n" 
		 "where m is the number of variants in the region", \
		 true, 0, "numeric");		
		TCLAP::ValueArg<std::string> priorFileArg("f","prior-file", \
			  "the file name specifying the prior probabilities of variants", \
								    true, "", "string");									
		TCLAP::SwitchArg outputStepwiseResultArg("s","stepwise", \
			"output stepwise result with rho confidence level", false);
		TCLAP::SwitchArg outputExhaustiveResultArg("e","exhaustive", \
	"output exhaustive search result with rho confidence level", false);
		TCLAP::SwitchArg outputExhaustiveStepwiseResultArg("x","mixed", \
        "output result first using exhaustive and then stepwise search "
		 "with rho confidence level", false);
		 
		cmd.xorAdd(priorArg, priorFileArg); 
		cmd.add(outputFileArg);
		cmd.add(outputExhaustiveStepwiseResultArg);
		cmd.add(outputExhaustiveResultArg);
		cmd.add(outputStepwiseResultArg);
		cmd.add(nSNPsArg);
		cmd.add(inputFileArg);
		 
		cmd.parse(argc, argv);

		inputFile = inputFileArg.getValue();
		outputFile = outputFileArg.getValue();
		nSNPs = nSNPsArg.getValue();
		if (priorArg.isSet()) {
			prior = priorArg.getValue();
			priorInAFile = false;
		} else if (priorFileArg.isSet()) {
			priorFile = priorFileArg.getValue();
			priorInAFile = true;
		} else {
			throw BadCommandLine();
		}
		outputStepwiseResult = outputStepwiseResultArg.getValue();
		outputExhaustiveResult = outputExhaustiveResultArg.getValue();
		outputExhaustiveStepwiseResult = \
		            outputExhaustiveStepwiseResultArg.getValue();

	} catch (TCLAP::ArgException &e)  {// catch any exceptions
		std::cerr << "error: " << e.error() << " for arg " \
		<< e.argId() << std::endl; 
	}
}

InputArgumentCAVIARBF::InputArgumentCAVIARBF() {
}

InputArgumentCAVIARBF::~InputArgumentCAVIARBF()
{
}

void InputArgumentCAVIARBF::processArguments(int argc, char** argv) {
	// Process the command line arguments to assign the following values
	//	  string zFilename;
	//	  string corFilename;
	//	  int priorType;
	//	  double priorValue;
	//	  int nSample;
	//	  int maxCausal;
	//	  string outputFile;

	try {  
		TCLAP::CmdLine cmd("Calculate Bayes factors based on "
							 "summary statistics", ' ', version);
		TCLAP::ValueArg<std::string> zFilenameArg("z","zfile", \
		                        "input file for marginal test statistics", \
								    true, "", "string");
		TCLAP::ValueArg<std::string> corFilenameArg("r","rfile", \
								"input file for correlation matrix", \
								true, "", "string");
		TCLAP::ValueArg<int> priorTypeArg("t","prior-type", \
			"prior type for variant effect size\n" 
			"\t0: specify sigmaa\n"
			"\t1: specify the proportion of variance explained (pve)", \
			true, 0, "integer");
		TCLAP::ValueArg<double> priorValueArg("a", "prior-value", \
							"the prior value associated with the prior type", \
								    true, 0, "numeric");	
		TCLAP::ValueArg<int> nSampleArg("n","sample-number", \
		                        "the total number of samples in the data", \
								    true, 0, "integer");	
		TCLAP::ValueArg<int> maxCausalArg("c","max-causal", \
						"the maximal number of causal variants in the model", \
								    true, 3, "integer");	
		TCLAP::ValueArg<std::string> outputFileArg("o","output", \
			  "the output file name for Bayes factors", \
								    true, "", "string");									
		 
		cmd.add(outputFileArg); 
		cmd.add(maxCausalArg);
		cmd.add(nSampleArg);
		cmd.add(priorValueArg);
		cmd.add(priorTypeArg);
		cmd.add(corFilenameArg);
		cmd.add(zFilenameArg);
		 
		cmd.parse(argc, argv);

		zFilename = zFilenameArg.getValue();
		corFilename = corFilenameArg.getValue();
		priorType = priorTypeArg.getValue();
		priorValue = priorValueArg.getValue();
		nSample = nSampleArg.getValue();
		maxCausal = maxCausalArg.getValue();
		outputFile = outputFileArg.getValue();

	} catch (TCLAP::ArgException &e)  {// catch any exceptions
		std::cerr << "error: " << e.error() << " for arg " \
		<< e.argId() << std::endl; 
	}
}