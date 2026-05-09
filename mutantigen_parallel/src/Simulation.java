/* Simulation functions, holds the host population */

import java.util.*;
import java.io.*;

//import com.javamex.classmexer.*;
import com.javamex.classmexer.*;

public class Simulation {

	// fields
	private List<HostPopulation> demes = new ArrayList<HostPopulation>();
	private double diversity;
	private double tmrca;
	private double netau;
	private double serialInterval;
	private double antigenicDiversity;
        private double meanLoad;
	
	private List<Double> diversityList = new ArrayList<Double>();
	private List<Double> tmrcaList = new ArrayList<Double>();	
	private List<Double> netauList = new ArrayList<Double>();		
	private List<Double> serialIntervalList = new ArrayList<Double>();		
	private List<Double> antigenicDiversityList = new ArrayList<Double>();
	private List<Double> nList = new ArrayList<Double>();
	private List<Double> sList = new ArrayList<Double>();	
	private List<Double> iList = new ArrayList<Double>();	
	private List<Double> rList = new ArrayList<Double>();		
	private List<Double> casesList = new ArrayList<Double>();
        
        private ArrayList<Double> daysList = new ArrayList<Double>();
        
        // Added these to track viral fitness dists
        private double meanRDist;
        private double varRDist;
        private double meanBetaDist;
        private double varBetaDist;
        private double meanSigmaDist;
        private double varSigmaDist;
        private double covBetaSigmaDist;
        
        private String simNum;
	
	// constructor
	public Simulation(String simNum) {
		this.simNum = simNum;
		for (int i = 0; i < Parameters.demeCount; i++) {
			if (Parameters.restartFromCheckpoint) {
				HostPopulation hp = new HostPopulation(i, true);
				demes.add(hp);
			}
			else {
				HostPopulation hp = new HostPopulation(i);
				demes.add(hp);
			}
		}
	}
	
	// methods
	
	public int getN() {
		int count = 0;
		for (int i = 0; i < Parameters.demeCount; i++) {
			HostPopulation hp = demes.get(i);
			count += hp.getN();
		}
		return count;
	}
	
	public int getS() {
		int count = 0;
		for (int i = 0; i < Parameters.demeCount; i++) {
			HostPopulation hp = demes.get(i);
			count += hp.getS();
		}
		return count;
	}	
	
	public int getI() {
		int count = 0;
		for (int i = 0; i < Parameters.demeCount; i++) {
			HostPopulation hp = demes.get(i);
			count += hp.getI();
		}
		return count;
	}	
	
	public int getR() {
		int count = 0;
		for (int i = 0; i < Parameters.demeCount; i++) {
			HostPopulation hp = demes.get(i);
			count += hp.getR();
		}
		return count;
	}		
	
	public int getCases() {
		int count = 0;
		for (int i = 0; i < Parameters.demeCount; i++) {
			HostPopulation hp = demes.get(i);
			count += hp.getCases();
		}
		return count;
	}	
	
	public double getDiversity() {
		return diversity;
	}		
	
	public double getNetau() {
		return netau;
	}	
	
	public double getTmrca() {
		return tmrca;
	}	
	
	public double getSerialInterval() {
		return serialInterval;	
	}		
	
	public double getAntigenicDiversity() {
		return antigenicDiversity;
	}		
		
	// proportional to infecteds in each deme
	public int getRandomDeme() {
		int n = Random.nextInt(0,getN()-1);
		int d = 0;
		int target = (demes.get(0)).getN();
		while (n < target) {
			d += 1;
			target += (demes.get(d)).getN();
		}
		return d;
	}
	
	// return random virus proportional to worldwide prevalence
	public Virus getRandomInfection() {
	
		Virus v = null;
		
		if (getI() > 0) {
	
			// get deme proportional to prevalence
			int n = Random.nextInt(0,getI()-1);
			int d = 0;
			int target = (demes.get(0)).getI();
			while (d < Parameters.demeCount) {
				if (n < target) {
					break;
				} else {
					d++;
					target += (demes.get(d)).getI();
				}	
			}
			HostPopulation hp = demes.get(d);
					
			// return random infection from this deme
			if (hp.getI()>0) {
				Host h = hp.getRandomHostI();
				v = h.getInfection();
			}
		
		}
		
		return v;
		
	}
        
        public double getMeanLoad() {
            HostPopulation hp = demes.get(0);
            double mean = hp.getMeanLoad();
            return mean;
        }
        
        public int getAnitgenicTypes() {
            HostPopulation hp = demes.get(0);
            int types = hp.getAntigenicCount();
            return types;
        }
        
        public int getCumlAntigenicTypes() {
            //int cumlTypes = AntigenicTree.nodeCount();
            ArrayList<Integer> treeTypes = AntigenicTree.getTreeTypes();
            int cumlTypes = treeTypes.size();
            return cumlTypes;
        }
        
//        public void reconstructAntDynamics(ArrayList<Integer> types) {
//            
//            //Print to file
//        }
        
	
	// return random host from random deme
	public Host getRandomHost() {
		int d = Random.nextInt(0,Parameters.demeCount-1);
		HostPopulation hp = demes.get(d);
		return hp.getRandomHost();
	}
	
	public double getAverageRisk(Phenotype p) {
		
		double averageRisk = 0;
		for (int i = 0; i < 10000; i++) {
			Host h = getRandomHost();
			Phenotype[] history = h.getHistory();
			averageRisk += p.riskOfInfection(history);
		}
		averageRisk /= 10000.0;
		return averageRisk;
		
	}
		
	public void printImmunity() {
	
		try {
			String newFileName = "outfiles/out_" + Mutantigen.globalSimNum + ".immunity";
			File immunityFile = new File(newFileName);
			immunityFile.delete();
			immunityFile.createNewFile();
			PrintStream immunityStream = new PrintStream(immunityFile);
			
			for (double x = VirusTree.xMin; x <= VirusTree.xMax; x += 0.5) {
				for (double y = VirusTree.yMin; y <= VirusTree.yMax; y += 0.5) {
				
					Phenotype p = PhenotypeFactory.makeArbitaryPhenotype(x,y);
					double risk = getAverageRisk(p);
					immunityStream.printf("%.4f,", risk);
				
				}
				immunityStream.println();
			}
			
			immunityStream.close();
		} catch(IOException ex) {
			System.out.println("Could not write to file"); 
			System.exit(0);
		}
	
	}
	
	public void printHostPopulation() {
	
		try {
			String newFileName = "outfiles/out_" + simNum + ".hosts";
			File hostFile = new File(newFileName);
			hostFile.delete();
			hostFile.createNewFile();
			PrintStream hostStream = new PrintStream(hostFile);
			for (int i = 0; i < Parameters.demeCount; i++) {
				HostPopulation hp = demes.get(i);
				hp.printHostPopulation(hostStream);
			}
			hostStream.close();
		} catch(IOException ex) {
			System.out.println("Could not write to file"); 
			System.exit(0);
		}
	
	}	
		
	public void makeTrunk() {
		for (int i = 0; i < Parameters.demeCount; i++) {
			HostPopulation hp = demes.get(i);
			hp.makeTrunk();
		}
	}	
	
	public void printState() {
                
		System.out.printf("%d\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%d\t%d\t%d\t%d\t%d\t%.3f\t%d\t%d\n", Parameters.day, getDiversity(), getTmrca(),  getNetau(), getSerialInterval(), getAntigenicDiversity(), getN(), getS(), getI(), getR(), getCases(), getMeanLoad(), getAnitgenicTypes(), getCumlAntigenicTypes());
		
		if (Parameters.memoryProfiling && Parameters.day % 100 == 0) {
			long noBytes = 0;
                        //long noBytes = MemoryUtil.deepMemoryUsageOf(this);
			//System.out.println("Total: " + noBytes);
			//HostPopulation hp = demes.get(1);
			//noBytes = MemoryUtil.deepMemoryUsageOf(hp);
			//System.out.println("One host population: " + noBytes);
			//HostPopulation hp = demes.get(0);
                        //Host h = hp.getRandomHostS();
			//long noBytes = MemoryUtil.deepMemoryUsageOf(h);
			//System.out.println("One susceptible host with " +  h.getHistoryLength() + " previous infection: " + noBytes);
			//h.printHistory();
			if (getI() > 0) {
				//Virus v = getRandomInfection();
				//noBytes = MemoryUtil.memoryUsageOf(v);
				//System.out.println("One virus: " + noBytes);
				noBytes = MemoryUtil.deepMemoryUsageOf(VirusTree.getTips());
                                double doubleBytes = (double) noBytes;
                                double megaBytes = 6000.0;
                                double conversion = 1050000.0;
                                double percent = doubleBytes / (megaBytes * conversion);
				System.out.println("Virus tree: " + noBytes + " " + "(" + percent + ")");
                                
                                noBytes = MemoryUtil.deepMemoryUsageOf(AntigenicTree.getDistArray());
                                doubleBytes = (double) noBytes;
                                percent = doubleBytes / (megaBytes * conversion);
				System.out.println("Antigenic tree: " + noBytes + " " + "(" + percent + ")");
			}
		}
		
	}

	public void printHeader(PrintStream stream) {
		stream.print("date\tdiversity\ttmrca\tnetau\tserialInterval\tantigenicDiversity\ttotalN\ttotalS\ttotalI\ttotalR\ttotalCases");
		for (int i = 0; i < Parameters.demeCount; i++) {
			HostPopulation hp = demes.get(i);
			hp.printHeader(stream);
		}
		stream.println();
	}
	
	public void printState(PrintStream stream) {
		stream.printf("%.4f\t%.4f\t%.4f\t%.4f\t%.5f\t%.4f\t%d\t%d\t%d\t%d\t%d", Parameters.getDate(), getDiversity(), getTmrca(), getNetau(), getSerialInterval(), getAntigenicDiversity(), getN(), getS(), getI(), getR(), getCases());
		for (int i = 0; i < Parameters.demeCount; i++) {
			HostPopulation hp = demes.get(i);
			hp.printState(stream);
		}
		stream.println();
	}
        
        public void printMutations(PrintStream stream) {
		//stream.printf("%.4f\t%.4f\t%.4f\t%.4f\t%.5f\t%.4f\t%d\t%d\t%d\t%d\t%d", Parameters.getDate(), getDiversity(), getTmrca(), getNetau(), getSerialInterval(), getAntigenicDiversity(), getN(), getS(), getI(), getR(), getCases());
		for (int i = 0; i < Parameters.demeCount; i++) {
			HostPopulation hp = demes.get(i);
			hp.printMutations(stream);
		}
		stream.println();
	}
        
        public void printFitness(PrintStream stream) {
		stream.printf("%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f", Parameters.getDate(), meanRDist, varRDist, meanBetaDist, varBetaDist, meanSigmaDist, varSigmaDist, covBetaSigmaDist);
		stream.println();
	}
        
        public void printFitnessSamples() {
            
                try {
                    String fileName = "outfiles/out" + simNum + ".fitnessSamples_t" + Double.toString(Parameters.day);
                    File samplesFile = new File(fileName);
                    samplesFile.delete();
                    samplesFile.createNewFile();
                    PrintStream samplesStream = new PrintStream(samplesFile);
                    
                    samplesStream.print("betas\tS/N\tR");
                    
                    // Get samples
                    HostPopulation hp = demes.get(0);
                    ArrayList<ArrayList<Double>> samplesData = hp.getViralFitnessSamples();
                    
                    int samples = samplesData.get(0).size();
                    for (int s = 0; s < samples; s++) {
                        double beta = samplesData.get(0).get(s);
                        double sigma = samplesData.get(1).get(s);
                        double R = samplesData.get(2).get(s);
                        samplesStream.printf("%.6f\t%.6f\t%.6f", beta, sigma, R);
                        samplesStream.println();
                    }
                    
                    samplesStream.close();
                    
                } catch(IOException ex) {
                    System.out.println("Could not write to file"); 
                    System.exit(0);
		}
            
        }
        
        public void printFitnessHeader(PrintStream stream) {
		stream.print("date\tmeanR\tvarR\tmeanBeta\tvarBeta\tmeanSigma\tvarSigma\tcovBetaSigma");
		stream.println();
	}
        
        public void printAntigenicDistances(ArrayList<Integer> types) {
                AntigenicTree.printDistanceMatrix(types);
        }
        
        public void printAntigenicDistancesBetweenTips() {
                AntigenicTree.printDistanceMatrixBetweenTips();
        }
	
	public void printSummary() {
	
		try {
			String newFileName = "outfiles/out_" + simNum + ".summary";
			File summaryFile = new File(newFileName);
			summaryFile.delete();
			summaryFile.createNewFile();
			PrintStream summaryStream = new PrintStream(summaryFile);
			summaryStream.printf("parameter\tfull\n");
			
			summaryStream.printf("diversity\t%.4f\n", mean(diversityList));
			summaryStream.printf("tmrca\t%.4f\n", mean(tmrcaList));
			summaryStream.printf("netau\t%.4f\n", mean(netauList));		
			summaryStream.printf("serialInterval\t%.5f\n", mean(serialIntervalList));	
			summaryStream.printf("antigenicDiversity\t%.4f\n", mean(antigenicDiversityList));	
			summaryStream.printf("N\t%.4f\n", mean(nList));		
			summaryStream.printf("S\t%.4f\n", mean(sList));		
			summaryStream.printf("I\t%.4f\n", mean(iList));		
			summaryStream.printf("R\t%.4f\n", mean(rList));		
			summaryStream.printf("cases\t%.4f\n", mean(casesList));					
			
			summaryStream.close();
		} catch(IOException ex) {
			System.out.println("Could not write to file"); 
			System.exit(0);
		}
	
	}
	
	private double mean(List<Double> list) {
		double mean = 0;
		if(!list.isEmpty()) {
			for (Double item : list) {
				mean += (double) item;
    		}
    	mean /= (double) list.size();
  		}
  		return mean;
	}	
		
	public void updateDiversity() {

		diversity = 0.0;
		tmrca = 0.0;
		antigenicDiversity = 0.0;		
		netau = 0.0;
		serialInterval = 0.0;
		
		double coalCount = 0.0;	
		double coalOpp = 0.0;
		double coalWindow = Parameters.netauWindow / 365.0;
		int sampleCount = Parameters.diversitySamplingCount;
		
		for (int i = 0; i < sampleCount; i++) {
			Virus vA = getRandomInfection();
			Virus vB = getRandomInfection();
			if (vA != null && vB != null) {
				double dist = vA.distance(vB);
				diversity += dist;
				if (dist > tmrca) {
					tmrca = dist;
				}
				antigenicDiversity += vA.antigenicDistance(vB);
				coalOpp += coalWindow;
				coalCount += vA.coalescence(vB, coalWindow);
				serialInterval += vA.serialInterval();
			}
		}	
	
		diversity /= (double) sampleCount;
		tmrca /= 2.0;	
		netau = coalOpp / coalCount;
		serialInterval /= (double) sampleCount;
		antigenicDiversity /= (double) sampleCount;			
		
	}
        
        public void updateFitnessDists() {
            
                HostPopulation hp = demes.get(0);
                //double mean = hp.getAntigenicLoadDistribution(); // is this used???
                VirusFitnessDist fitDist = hp.getViralFitnessDistribution(); // VirusFitnessDist holds variables
                meanRDist = fitDist.meanR;
                varRDist = fitDist.varR;
                meanBetaDist = fitDist.meanBeta;
                varBetaDist = fitDist.varBeta;
                meanSigmaDist = fitDist.meanSigma;
                varSigmaDist = fitDist.varSigma;
                covBetaSigmaDist = fitDist.covBetaSigma;
                
                
        }
	
	public void pushLists() {
		diversityList.add(diversity);
		tmrcaList.add(tmrca);
		netauList.add(netau);
		serialIntervalList.add(serialInterval);
		antigenicDiversityList.add(antigenicDiversity);
		nList.add((double) getN());
		sList.add((double) getS());
		iList.add((double) getI());
		rList.add((double) getR());
		casesList.add((double) getCases());		
	}
				
	public void resetCases() {
		for (int i = 0; i < Parameters.demeCount; i++) {	
			HostPopulation hp = demes.get(i);
			hp.resetCases();
		}
	}
		
	public void stepForward() {
				
		for (int i = 0; i < Parameters.demeCount; i++) {		
			HostPopulation hp = demes.get(i);
			
                        //hp.stepForward(); // Trevor's method
                        
                        hp.stepForwardMutAtTrans(); // My method with mutationa at transmission
			
                        for (int j = 0; j < Parameters.demeCount; j++) {
				if (i != j) {
					HostPopulation hpOther = demes.get(j);
					hp.betweenDemeContact(hpOther);
				}
			}
		}
							
		Parameters.day++;
		
	}
	
	public void run() {
	
		try {
			String newFileName = "outfiles/out_" + simNum + ".timeseries";
			File seriesFile = new File(newFileName);		
			seriesFile.delete();
			seriesFile.createNewFile();
			PrintStream seriesStream = new PrintStream(seriesFile);
			System.out.println("day\tdiversity\ttmrca\tnetau\tserialInterval\tantigenicDiversity\tN\tS\tI\tR\tcases\tmeanLoad\tantigenicTypes\tcumulativeTypes");
			printHeader(seriesStream);
                        
                        String newNewFileName = "outfiles/out_" + simNum + ".mutationSeries";
                        File mutationFile = new File(newNewFileName);
                        mutationFile.delete();
                        mutationFile.createNewFile();
                        PrintStream mutationStream = new PrintStream(mutationFile);
			
		       String lastFileName = "outfiles/out_" + simNum + ".viralFitnessSeries";
                        File fitnessFile = new File(lastFileName);
                        fitnessFile.delete();
                        fitnessFile.createNewFile();
                        PrintStream fitnessStream = new PrintStream(fitnessFile);
                        printFitnessHeader(fitnessStream);
                        
                        
			for (int i = 0; i < Parameters.endDay; i++) {
				
				stepForward();
				
				if (Parameters.day % Parameters.printStep == 0) {
                                    daysList.add(Parameters.getDate());
                                    updateDiversity();
                                    updateFitnessDists(); // new for tracking viral fitness distributions
                                    printState();
                                    if (Parameters.day > Parameters.burnin) {
                                        printState(seriesStream);
                                        printMutations(mutationStream);
                                        printFitness(fitnessStream);
                                        pushLists();
                                    }
                                    resetCases();
				}
                                
                                if (Parameters.day % Parameters.printFitSamplesStep == 0) {
                                    // Print list of randomly sampled virus's beta, S/N and R values
                                    printFitnessSamples();
                                }
				
				if (getI()==0) {
					if (Parameters.repeatSim) {
						reset();
						i = 0; 
						seriesFile.delete();
						seriesFile.createNewFile();
						seriesStream = new PrintStream(seriesFile);
						printHeader(seriesStream);
					} else {
						break;
					}
				}
			}
			
			seriesStream.close();
                        mutationStream.close();
                        fitnessStream.close();
		} catch(IOException ex) {
			System.out.println("Could not write to file"); 
			System.exit(0);
		}	
		
		// Summary
		printSummary();
					
		// tree reduction
		
		double countsTemp = getR();
		int paramsDemeCount = Parameters.demeCount;
		// int numOfTips = VirusTree.getTips().size();
		int numOfTips = 0;
		
			for (Virus v : VirusTree.getTips()) {
                                        if ((v.getBirth()*365.25) >= Parameters.tipSamplingStartDay & (v.getBirth()*365.25) <= Parameters.tipSamplingEndDay) {
                                           numOfTips+=1;
                                        }
			}
		
		System.out.printf("The total number of tips reported in VirusTree between sampling days: %d \n", numOfTips);
		//double demeCountTemp = (double) VirusTree.getDemeCount(d);
		//double keepProportion = (double) Parameters.tipSamplesPerDeme / (double) VirusTree.getDemeCount(d);
		
		VirusTree.pruneTips();
		VirusTree.markTips();		
	
		// tree prep
		makeTrunk();
		VirusTree.fillBackward();                                       // assign parents children			
		VirusTree.sortChildrenByDescendants();                          // sort children by number of descendents
		VirusTree.setLayoutByDescendants();
		
		VirusTree.getMRCASeries(daysList);                              // prints out time of MRCA in tree over time
		
		String filePath5 = "logfiles/manual_output_" + Mutantigen.globalSimNum + ".txt";
        	String textToAppend5 = "About to streamline the virus tree.\n";

        try (FileWriter fw = new FileWriter(filePath5, true)) {
            fw.write(textToAppend5);
        } catch (IOException e) {
            System.err.println("An error occurred: " + e.getMessage());
        }  
        
        
                VirusTree.streamline();                                         // collapses nodes with single descentdent in tree			
		
                
                // reconstruct dynamics of antigenic types in tree
             VirusTree.getTreeTypes();   
             String filePath = "logfiles/manual_output_" + Mutantigen.globalSimNum + ".txt";
        	String textToAppend = "About to reconstruct ant dynamics.\n";

            try (FileWriter fw = new FileWriter(filePath, true)) {
                fw.write(textToAppend);
            } catch (IOException e) {
                System.err.println("An error occurred: " + e.getMessage());
            }                                    // get all antigenic types in tree (including internals)
                //VirusTree.reconstructAntDynamics(daysList);                     // reconstruct the population dynamics of all antigenic types in tree
                //printAntigenicDistances(VirusTree.treeTypes);                   // print antigenic distance matrix for antigenic types in tree
                
    		// rotation
    		if (Parameters.pcaSamples) {
    			VirusTree.rotate();
    			VirusTree.flip();
    		}
        	 String filePath6 = "logfiles/manual_output_" + Mutantigen.globalSimNum + ".txt";
            String textToAppend6 = "About to write a bunch of stuff about the virus tree.\n";
    
            try (FileWriter fw = new FileWriter(filePath6, true)) {
                fw.write(textToAppend6);
            } catch (IOException e) {
                System.err.println("An error occurred: " + e.getMessage());
            }
        
        
		// tip and tree output
		VirusTree.printTips();			
		VirusTree.printBranches();	
                VirusTree.printNewick();                                        // Print basic newick tree
                            String filePath2 = "logfiles/manual_output_" + Mutantigen.globalSimNum + ".txt";
        	String textToAppend2 = "About to save the print sim map load file.\n";

        try (FileWriter fw = new FileWriter(filePath2, true)) {
            fw.write(textToAppend2);
        } catch (IOException e) {
            System.err.println("An error occurred: " + e.getMessage());
        }
        
                VirusTree.printSimMapLoad();                                    // Print SimMap tree with annotation for mutation load along lineages
                
                            String filePath3 = "logfiles/manual_output_" + Mutantigen.globalSimNum + ".txt";
        	String textToAppend3 = "About to save the call the method that prints the sim map antigenic.\n";

        try (FileWriter fw = new FileWriter(filePath3, true)) {
            fw.write(textToAppend3);
        } catch (IOException e) {
            System.err.println("An error occurred: " + e.getMessage());
        }
                
                VirusTree.printSimMapAntigenic();                               // Print SimMap tree with annotation for antigenic type along lineages
		
                VirusTree.printTrunkAntigenicShifts();                          // Print size of all antigenic transitions along trunk
                
                
                // Trevor's stuff:
                
		// mk output
		VirusTree.printMK();
		
		// immunity output
		if (Parameters.phenotypeSpace == "geometric") {
			VirusTree.updateRange();
			VirusTree.printRange();
			if (Parameters.immunityReconstruction) {
				printImmunity();
			}
		}		
		
		// detailed output
		if (Parameters.detailedOutput) {
			printHostPopulation();
		}	
		
	}
	
	public void reset() {
		Parameters.day = 0;
		diversity = 0;
		for (int i = 0; i < Parameters.demeCount; i++) {
			HostPopulation hp = demes.get(i);
			hp.reset();
		}
		VirusTree.clear();
	}

}
