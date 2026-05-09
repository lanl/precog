/* Stores parameters for use across simulation */
/* Start with parameters in source, implement input file later */
/* A completely static class.  */

import java.util.*;
import static java.lang.Math.*;
import java.io.*;

public class Parameters {
	
	// global parameters
	public static int day = 0;
	public static Virus urVirus = null;
	public static Phenotype urImmunity = null;
        public static double meanLoad = 0.0;
	
	// simulation parameters
	public static int burnin = 0;
	public static int endDay = 5000; 
	public static int printStep = 10;									// print to out.timeseries every week
	public static int tipSamplingStartDay = 0;
        public static int tipSamplingEndDay = 0;
        public static double tipSamplingRate = 0.0002;						// in samples per deme per day
	public static int tipSamplesPerDeme = 1000;
	public static boolean tipSamplingProportional = true;				// whether to sample proportional to prevalance
	public static double treeProportion = 0.1;							// proportion of tips to use in tree reconstruction
	public static int diversitySamplingCount = 1000;					// how many samples to draw to calculate diversity, Ne*tau, serial interval
	public static int netauWindow = 100;								// window in days to calculate Ne*tau		
	public static boolean repeatSim = true;								// repeat simulation until endDay is reached?
	public static boolean immunityReconstruction = false;				// whether to print immunity reconstruction to out.immunity
	public static boolean memoryProfiling = false;						// requires -javaagent:classmexer.jar to run
	public static double yearsFromMK = 1.0;
	public static boolean pcaSamples = false;							// whether to rotate and flip virus tree
	public static boolean detailedOutput = false;						// whether to output out.hosts and out.viruses files enabling checkpointing
	public static boolean restartFromCheckpoint = false;				// whether to load population from out.hosts

	// metapopulation parameters
	public static int demeCount = 3;
	public static String[] demeNames = {"north", "tropics", "south"};
	public static int[] initialNs = {1000000,1000000,1000000};	
	
	// host parameters
	public static double birthRate = 0.000091;				// in births per individual per day, 1/30 years = 0.000091
	public static double deathRate = 0.000091;				// in deaths per individual per day, 1/30 years = 0.000091
	public static boolean swapDemography = true;				// whether to keep overall population size constant
		
	// epidemiological parameters
	public static int initialI = 10;							// in individuals
	public static int initialDeme = 2;						// index of deme where infection starts, 1..n
	public static double initialPrR = 0.5; 					// as proportion of population
	public static double beta = 0.36; // 0.3					// in contacts per individual per day
	public static double nu = 0.2; //0.2						// in recoveries per individual per day
	public static double betweenDemePro = 0.0005;				// relative to within-deme beta	
        public static double externalMigration = 0.0;                           // relative to within-deme beta
        
	// transcendental immunity
	public static boolean transcendental = false;
	public static double immunityLoss = 0.01;					// in R->S per individual per day
	public static double initialPrT = 0.1;
        
        public static boolean backgroundImmunity = false;
	public static double backgroundDistance = 0.0;
        
	// seasonal betas
	public static double[] demeBaselines = {1,1,1};
	public static double[] demeAmplitudes = {0.1,0,0.1};
	public static double[] demeOffsets = {0,0,0.5};				// relative to the year
	
	// phenotype parameters
	public static String phenotypeSpace = "geometric";			// options include: "geometric", "geometric3d", "geometric10d"
	public static double muPhenotype = 0.005; 					// in mutations per individual per day

	// parameters specific to GeometricPhenotype
	public static double smithConversion = 0.1;					// multiplier to distance to give cross-immunity
	public static double homologousImmunity = 0.05;				// immunity raised to antigenically identical virus
	public static double initialTraitA = -6;	
	public static double meanStep = 0.3; 
	public static double sdStep = 0.3; 
	public static boolean mut2D = false;                                    // whether to mutate in a full 360 degree arc
		
	// parameters specific to MutLoadPhenotype
        public static double lambda;                                            // deleterious mutation rate
        public static double mutCost;                                           // deleterious mutation fitness cost
        public static double probLethal;                                        // probability of a lethal mutation
        public static double epsilon;                                           // prob beneficial mutation
        public static double epsilonSlope;                                      // dependece of epsilon on mean load
        
        // parameters specific to additive antigenic evolution
        public static double lambdaAntigenic;                                   // antigenic mutation rate
        public static double meanAntigenicSize;                                 // mean size of antigenic distributions
        public static double antigenicGammaShape;                               // shape parameter for antigenic mutations
        public static double thresholdAntigenicSize;
        public static int antigenicEvoStartDay;
        public static double cleanUpDistance;
        
        public static double demoNoiseScaler;
        
        // parameter for sampling host immunie histories
        public static int hostImmuneHistorySampleCount;
        public static int fitSampleCount;
        public static int printFitSamplesStep;
        
        // measured in years, starting at burnin
	public static double getDate() {
		return ((double) day - (double) burnin ) / 365.0;
	}
	
	public static double getSeasonality(int index) {
		double baseline = demeBaselines[index];
		double amplitude = demeAmplitudes[index];
		double offset = demeOffsets[index];
		double beta = baseline + amplitude * Math.cos(2*Math.PI*getDate() + 2*Math.PI*offset);
		return beta;
	}
		
	// initialize
	public static void initialize() {
            
		urVirus = new Virus();
		urImmunity = PhenotypeFactory.makeHostPhenotype();
                //System.out.println();
	}
	
	// load parameters.yml	
	public static void load(Object inputStr) {
		
		try {
			String strName = inputStr.toString();
			
			org.yaml.snakeyaml.Yaml yaml = new org.yaml.snakeyaml.Yaml();
			Map map = null;			
			InputStream input = new FileInputStream(new File(strName));
			map = (Map) yaml.load(input);
			input.close();
			
			System.out.println("Loading parameters from parameters.yml");
 
                        
                        burnin = toInteger(map.get("burnin")); //burnin = (int) map.get("burnin");
			endDay = toInteger(map.get("endDay")); //endDay = (int) map.get("endDay");
			printStep = toInteger(map.get("printStep")); //printStep = (int) map.get("printStep");
			tipSamplingStartDay = toInteger(map.get("tipSamplingStartDay"));
                        tipSamplingEndDay = toInteger(map.get("tipSamplingEndDay"));
                        tipSamplingRate = toDouble(map.get("tipSamplingRate")); //tipSamplingRate = (double) map.get("tipSamplingRate");
			tipSamplesPerDeme = toInteger(map.get("tipSamplesPerDeme")); //tipSamplesPerDeme = (int) map.get("tipSamplesPerDeme");
			tipSamplingProportional = toBoolean(map.get("tipSamplingProportional")); //tipSamplingProportional = (boolean) map.get("tipSamplingProportional");
			treeProportion = toDouble(map.get("treeProportion")); //treeProportion = (double) map.get("treeProportion");
			diversitySamplingCount = toInteger(map.get("diversitySamplingCount")); //diversitySamplingCount = (int) map.get("diversitySamplingCount");	
			netauWindow = toInteger(map.get("netauWindow")); //netauWindow = (int) map.get("netauWindow");	
			repeatSim = toBoolean(map.get("repeatSim")); //repeatSim = (boolean) map.get("repeatSim");
			immunityReconstruction = toBoolean(map.get("immunityReconstruction")); //immunityReconstruction = (boolean) map.get("immunityReconstruction");
			memoryProfiling = toBoolean(map.get("memoryProfiling")); //memoryProfiling = (boolean) map.get("memoryProfiling");
			yearsFromMK = toDouble(map.get("yearsFromMK")); //yearsFromMK = (double) map.get("yearsFromMK");
			pcaSamples = toBoolean(map.get("pcaSamples")); //pcaSamples = (boolean) map.get("pcaSamples");
			detailedOutput = toBoolean(map.get("detailedOutput")); //detailedOutput = (boolean) map.get("detailedOutput");
			restartFromCheckpoint = toBoolean(map.get("restartFromCheckpoint")); //restartFromCheckpoint = (boolean) map.get("restartFromCheckpoint");
			demeCount = toInteger(map.get("demeCount")); //demeCount = (int) map.get("demeCount");
			demeNames = toStringArray((List<String>) map.get("demeNames"));
			initialNs = toIntArray((List<Integer>) map.get("initialNs"));		
			birthRate = toDouble(map.get("birthRate")); //birthRate = (double) map.get("birthRate");
			deathRate = toDouble(map.get("deathRate")); //deathRate = (double) map.get("deathRate");
			swapDemography = toBoolean(map.get("swapDemography")); //swapDemography = (boolean) map.get("swapDemography");
			initialI = toInteger(map.get("initialI")); //initialI = (int) map.get("initialI");
			initialDeme = toInteger(map.get("initialDeme")); //initialDeme = (int) map.get("initialDeme");
			initialPrR = toDouble(map.get("initialPrR")); //initialPrR = (double) map.get("initialPrR");
			beta = toDouble(map.get("beta")); //beta = (double) map.get("beta");
			nu = toDouble(map.get("nu")); //nu = (double) map.get("nu");
			betweenDemePro = toDouble(map.get("betweenDemePro")); //betweenDemePro = (double) map.get("betweenDemePro");
			externalMigration = toDouble(map.get("externalMigration"));
                        transcendental = toBoolean(map.get("transcendental")); //transcendental = (boolean) map.get("transcendental");
			immunityLoss = toDouble(map.get("immunityLoss")); //immunityLoss = (double) map.get("immunityLoss");
			initialPrT = toDouble(map.get("initialPrT")); //initialPrT = (double) map.get("initialPrT");
			backgroundImmunity = toBoolean(map.get("backgroundImmunity"));
                        backgroundDistance = toDouble(map.get("backgroundDistance"));
                        demeBaselines = toDoubleArray((List<Double>) map.get("demeBaselines"));	
			demeAmplitudes = toDoubleArray((List<Double>) map.get("demeAmplitudes"));
			demeOffsets = toDoubleArray((List<Double>) map.get("demeOffsets"));		
			phenotypeSpace = (String) map.get("phenotypeSpace");
			muPhenotype = toDouble(map.get("muPhenotype")); //muPhenotype = (double) map.get("muPhenotype");
			smithConversion = toDouble(map.get("smithConversion")); //smithConversion = (double) map.get("smithConversion");
			homologousImmunity = toDouble(map.get("homologousImmunity")); //homologousImmunity = (double) map.get("homologousImmunity");
			initialTraitA = toDouble(map.get("initialTraitA")); //initialTraitA = (double) map.get("initialTraitA");
			meanStep = toDouble(map.get("meanStep")); //meanStep = (double) map.get("meanStep");
			sdStep = toDouble(map.get("sdStep")); //sdStep = (double) map.get("sdStep");
			mut2D = toBoolean(map.get("mut2D")); //mut2D = (boolean) map.get("mut2D");
                        
                        // added these for mutational load
                        lambda = toDouble(map.get("lambda"));
                        mutCost = toDouble(map.get("mutCost"));
                        probLethal = toDouble(map.get("probLethal"));
                        epsilon = toDouble(map.get("epsilon"));
                        epsilonSlope = toDouble(map.get("epsilonSlope"));
                        
                        // added this for additive antigenic evolution
                        lambdaAntigenic = toDouble(map.get("lambdaAntigenic"));
                        meanAntigenicSize = toDouble(map.get("meanAntigenicSize"));
                        antigenicGammaShape = toDouble(map.get("antigenicGammaShape"));
                        thresholdAntigenicSize = toDouble(map.get("thresholdAntigenicSize"));
                        antigenicEvoStartDay = toInteger(map.get("antigenicEvoStartDay"));
                        cleanUpDistance = toDouble(map.get("cleanUpDistance"));
                        
                        demoNoiseScaler = toDouble(map.get("demoNoiseScaler"));
                        
                        hostImmuneHistorySampleCount = toInteger(map.get("hostImmuneHistorySampleCount"));
                        
                        fitSampleCount = toInteger(map.get("fitSampleCount"));
                        printFitSamplesStep = toInteger(map.get("printFitSamplesStep"));

		
		} catch (IOException e) {
			System.out.println("Cannot load parameters.yml, using defaults");
		}		
	
	}
                
	//I added this one
        private static int toInteger(Object input) {
                String inStr = input.toString();
                int ret = Integer.valueOf(inStr);
                return ret;    
        }
        
        //I added this one
        private static double toDouble(Object input) {
                String inStr = input.toString();
                double ret = Double.valueOf(inStr);
                return ret;
        }
        
        //I added this one
        private static boolean toBoolean(Object input) {
                String inStr = input.toString();
                boolean ret = Boolean.valueOf(inStr);
                return ret;
        }
        
        
	private static int[] toIntArray(List<Integer> list){
  		int[] ret = new int[list.size()];
  		for (int i = 0; i < ret.length; i++) {
    		ret[i] = list.get(i);
    	}
  		return ret;
	}
	
	private static double[] toDoubleArray(List<Double> list){
  		double[] ret = new double[list.size()];
  		for (int i = 0; i < ret.length; i++) {
    		ret[i] = list.get(i);
    	}
  		return ret;
	}	
	
	private static String[] toStringArray(List<String> list){
  		String[] ret = new String[list.size()];
  		for (int i = 0; i < ret.length; i++) {
    		ret[i] = list.get(i);
    	}
  		return ret;
	}	
	
}