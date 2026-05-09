/* A population of host individuals */

import java.util.*;
import java.io.*;
import java.util.regex.*;

public class HostPopulation {

	// fields
	private int deme;
	private String name;	
	private int cases;	
	private List<Host> susceptibles = new ArrayList<Host>();
	private List<Host> infecteds = new ArrayList<Host>();	
	private List<Host> recovereds = new ArrayList<Host>();		// this is the transcendental class, immune to all forms of virus  
	private double diversity;
	private double tmrca;
	private double netau;	
	private double serialInterval;
	private double antigenicDiversity;
        
        private int infectionsNow;
        private int recoveriesNow;
        private List<Host> newParents;
        private List<Host> newInfections;
        
        //private double meanLoad;
        
        ArrayList<ArrayList<Host>> infectionsByClass;

	// construct population, using Virus v as initial infection
	public HostPopulation(int d) {
	
		// basic parameters
		deme = d;
		name = Parameters.demeNames[deme];
		int initialR = 0;
		if (Parameters.transcendental) {
			initialR = (int) ((double) Parameters.initialNs[deme] * Parameters.initialPrT);
		}
	
		// fill population with susceptibles
		int initialS = Parameters.initialNs[deme] - initialR;
		if (deme == Parameters.initialDeme - 1) {
			initialS -= Parameters.initialI;
		}
		for (int i = 0; i < initialS; i++) {
			Host h = new Host();
			susceptibles.add(h);
		}
		
		// fill population with recovereds
		for (int i = 0; i < initialR; i++) {
			Host h = new Host();			
			recovereds.add(h);
		}		
		
                double theta = Parameters.lambda * (1-Parameters.probLethal) / Parameters.mutCost;
		if (deme == Parameters.initialDeme - 1) {
		
			// infect some individuals
			for (int i = 0; i < Parameters.initialI; i++) {
				
                                //Start with a homogenous population
                                Virus v = new Virus(Parameters.urVirus, deme);
                                Phenotype p = new MutLoadPhenotype(theta);
                                v.setPhenotype(p);
                                Host h = new Host(v);
				infecteds.add(h);
			}	
		
		}
		
	}
	
	// construct checkpointed host population and infecting viruses
	public HostPopulation(int d, boolean checkpoint) {
	
		if (checkpoint == true) {
		
			deme = d;
			name = Parameters.demeNames[deme];
		
			try {
    			BufferedReader in = new BufferedReader(new FileReader("out.hosts"));
    			String line;
    			while ((line = in.readLine()) != null) {
    				Pattern regex = Pattern.compile(":");
    				String[] items = regex.split(line);
    				int thisDeme = Integer.parseInt(items[0]);
    				String sVirus = items[1];
    				String sHist = items[2];
        			if (thisDeme == deme) {
        				Host h = new Host(deme, sVirus, sHist);
        				if (sVirus.equals("n")) {
        					susceptibles.add(h);	
        				}
        				else {
        					infecteds.add(h);
        				}
        			}
    			}
    		in.close();
			} 
			catch (IOException ex) {
				System.out.println("Could not read in out.hosts"); 
				System.exit(0);
			}
		
		}
	
	}
	
	// accessors
	public int getN() {
		return susceptibles.size() + infecteds.size() + recovereds.size();
	}
	public int getS() {
		return susceptibles.size();
	}
	public int getI() {
		return infecteds.size();
	}
	public int getR() {
		return recovereds.size();
	}	
	public double getPrS() {
		return (double) getS() / (double) getN();
	}
	public double getPrI() {
		return (double) getI() / (double) getN();
	}
	public double getPrR() {
		return (double) getR() / (double) getN();
	}	
	public int getRandomN() {
		return Random.nextInt(0,getN()-1);
	}
	public int getRandomS() {
		return Random.nextInt(0,getS()-1);
	}
	public int getRandomI() {
		return Random.nextInt(0,getI()-1);
	}
	public int getRandomR() {
		return Random.nextInt(0,getR()-1);
	}
	
	public Host getRandomHost() {
		// figure out whether to pull from S, I or R
		Host h = null;
		double n = Random.nextDouble(0.0,1.0);
		if (n < getPrS()) {
			h = getRandomHostS();
		}
		else if (n > getPrS() && n < getPrS() + getPrI()) {
			h = getRandomHostI();
		}
		else if (n > getPrS() + getPrI()) {
			h = getRandomHostR();
		}
		return h;
	}
	
	public Host getRandomHostS() {
		int index = Random.nextInt(0,getS()-1);
		return susceptibles.get(index);
	}
	public Host getRandomHostI() {
		Host h = null;
		if (getI() > 0) {
			int index = Random.nextInt(0,getI()-1);
			h = infecteds.get(index);
		}
		return h;
	}
	public Host getRandomHostR() {
		Host h = null;
		if (getR() > 0) {	
			int index = Random.nextInt(0,getR()-1);
			h = recovereds.get(index);
		}
		return h;
	}	
	
	public Virus getRandomInfection() {
		Virus v = null;
		Host h = getRandomHostI();
		if (h != null) {
			v = h.getInfection();
		}
		return v;
	}	
	
	public void resetCases() {
		cases = 0;
	}
	public int getCases() {
		return cases;
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
        
        public double getMeanLoad() {
            return Parameters.meanLoad;
        }
        
        public void computeMeanLoad() {
            
            double mean = 0;
            int infections = getI();
            if (infections > 0) {
                int totalLoad = 0;
                for (int i = 0; i < infections; i++) {
                    Host h = infecteds.get(i);
                    Virus v = h.getInfection();
                    Phenotype p = v.getPhenotype();
                    totalLoad += p.mutLoad();
                }
                mean = totalLoad / infections;
            }
            Parameters.meanLoad = mean;
            
            //return meanLoad;
            
        }
        
        public VirusFitnessDist getViralFitnessDistribution() {
            
            VirusFitnessDist fitDist = new VirusFitnessDist();
                
            // Get mean immunity in host population to each unique type
            ArrayList<Integer> typeList = getAntigenicTypes(); // returns list of all unique antigenic types
            ArrayList<Double> typeImmunities = getMeanTypeImmunities(typeList); //store average immunity in host population fotyper each 
            
            // Get beta, sigma and R distributions in viral population
            ArrayList<Double> betas = new ArrayList<Double>();
            ArrayList<Double> sigmas = new ArrayList<Double>();
            ArrayList<Double> Rs = new ArrayList<Double>();
            
            double currS = getS();
            double currN = getN();
            double currSOverN = currS / currN;
            //double constants = 1/(Parameters.birthRate + Parameters.nu);
            
            //final double constantRTerm = currSOverN * (1/(Parameters.birthRate + Parameters.nu));
            final double constantRTerm = 1/(Parameters.birthRate + Parameters.nu);
            int infections = getI();
            
            double sumBetas = 0;
            double sumSigmas = 0;
            double sumRs = 0;
            
            boolean logValues = true;
            
            for (int i = 0; i < infections; i++) {
                
                Host h = infecteds.get(i);
                Virus v = h.getInfection();
                Phenotype p = v.getPhenotype();
                
                // Betas
                final int k = p.mutLoad();
                double betaK = Parameters.beta * Math.pow((1-Parameters.mutCost), k);
                if (logValues) {
                    betaK = Math.log(betaK);
                }
                sumBetas += betaK;
                betas.add(betaK);
                
                // Sigmas
                int typeLoc = typeList.indexOf(p.antigenicType());
                double sigma = typeImmunities.get(typeLoc) * currSOverN; // multiplied by S/N to get effective S/N
                if (logValues) {
                    sigma = Math.log(sigma);
                }
                sumSigmas += sigma;
                sigmas.add(sigma);
                
                double viralR = 0.0;
                if (logValues) {
                    viralR = betaK + sigma + Math.log(constantRTerm);
                } else {
                    viralR = betaK * sigma * constantRTerm;
                }
                sumRs += viralR;
                Rs.add(viralR);
                
            }
            
            // Compute mean and variance of betas
            double meanBeta = sumBetas / infections;
            double sumSqrDevsBeta = 0.0;
            for (int i = 0; i < infections; i++) {
                final double sqrDev = Math.pow((betas.get(i) - meanBeta),2);
                sumSqrDevsBeta += sqrDev;
            }
            final double varBeta = sumSqrDevsBeta / infections;
            fitDist.meanBeta = meanBeta;
            fitDist.varBeta = varBeta;

            // Compute mean and variance of sigmas
            double meanSigma = sumSigmas / infections;
            double sumSqrDevsSigmas = 0.0;
            for (int i = 0; i < infections; i++) {
                final double sqrDev = Math.pow((sigmas.get(i) - meanSigma),2);
                sumSqrDevsSigmas += sqrDev;
            }
            final double varSigma = sumSqrDevsSigmas / infections;
            fitDist.meanSigma = meanSigma;
            fitDist.varSigma = varSigma;
            
            // Compute mean and variance of Rs
            double meanR = sumRs / infections;
            double sumSqrDevsR = 0.0;
            for (int i = 0; i < infections; i++) {
                final double sqrDev = Math.pow((Rs.get(i) - meanR),2);
                sumSqrDevsR += sqrDev;
            }
            final double varR = sumSqrDevsR / infections;
            fitDist.meanR = meanR;
            fitDist.varR = varR;
            
            // Compute covariance of beta and sigma
            double sumProdDevs = 0.0;
            for (int i = 0; i < infections; i++) {
                final double prodDev = (betas.get(i) - meanBeta) * (sigmas.get(i) - meanSigma);
                sumProdDevs += prodDev;
            }
            final double covBetaSigma = sumProdDevs / infections;
            fitDist.covBetaSigma = covBetaSigma;
            
            return fitDist;
            
            
        }
        
        public ArrayList<ArrayList<Double>> getViralFitnessSamples() {
            
            ArrayList<ArrayList<Double>> samplesData = new ArrayList<ArrayList<Double>>();
            
            ArrayList<Integer> typeList = getAntigenicTypes(); // returns list of all unique antigenic types
            ArrayList<Double> typeImmunities = getMeanTypeImmunities(typeList); //store average immunity in host population fotyper each 
            
            int infections = getI();
            int samples;
            if (infections > Parameters.fitSampleCount) {
                samples = infections; 
            } else {
                samples = Parameters.fitSampleCount;
            }
            
            double currS = getS();
            double currN = getN();
            double currSOverN = currS / currN;
            final double constantRTerm = 1/(Parameters.birthRate + Parameters.nu);
            
            // Get beta, sigma and R distributions in viral population
            ArrayList<Double> betas = new ArrayList<Double>();
            ArrayList<Double> sigmas = new ArrayList<Double>();
            ArrayList<Double> Rs = new ArrayList<Double>();
            
            for (int i = 0; i < samples; i++) {
                
                Host h = getRandomHostI(); //infecteds.get(i);
                Virus v = h.getInfection();
                Phenotype p = v.getPhenotype();
                
                // Betas
                final int k = p.mutLoad();
                double betaK = Parameters.beta * Math.pow((1-Parameters.mutCost), k);
                betas.add(betaK);
                
                // Sigmas
                int typeLoc = typeList.indexOf(p.antigenicType());
                double sigma = typeImmunities.get(typeLoc) * currSOverN;
                sigmas.add(sigma);
                
                // Rs
                double viralR = betaK * sigma * constantRTerm;
                Rs.add(viralR);
                
            }
            
            samplesData.add(betas);
            samplesData.add(sigmas);
            samplesData.add(Rs);
            
            return samplesData;
            
        }
        
        // Old version where we were tracking sigma instead of S/N 
        public VirusFitnessDist getViralFitnessDistributionOld() {
            
            VirusFitnessDist fitDist = new VirusFitnessDist();
                
            // Get mean immunity in host population to each unique type
            ArrayList<Integer> typeList = getAntigenicTypes(); // returns list of all unique antigenic types
            ArrayList<Double> typeImmunities = getMeanTypeImmunities(typeList); //store average immunity in host population fotyper each 
            
            // Get beta, sigma and R distributions in viral population
            ArrayList<Double> betas = new ArrayList<Double>();
            ArrayList<Double> sigmas = new ArrayList<Double>();
            ArrayList<Double> Rs = new ArrayList<Double>();
            
            double currS = getS();
            double currN = getN();
            double currSOverN = currS / currN;
            //double constants = 1/(Parameters.birthRate + Parameters.nu);
            
            final double constantRTerm = currSOverN * (1/(Parameters.birthRate + Parameters.nu));
            int infections = getI();
            
            double sumBetas = 0;
            double sumSigmas = 0;
            double sumRs = 0;
            
            boolean logValues = true;
            
            for (int i = 0; i < infections; i++) {
                
                Host h = infecteds.get(i);
                Virus v = h.getInfection();
                Phenotype p = v.getPhenotype();
                
                // Betas
                final int k = p.mutLoad();
                double betaK = Parameters.beta * Math.pow((1-Parameters.mutCost), k);
                if (logValues) {
                    betaK = Math.log(betaK);
                }
                sumBetas += betaK;
                betas.add(betaK);
                
                // Sigmas
                int typeLoc = typeList.indexOf(p.antigenicType());
                double sigma = typeImmunities.get(typeLoc);
                if (logValues) {
                    sigma = Math.log(sigma);
                }
                sumSigmas += sigma;
                sigmas.add(sigma);
                
                double viralR = 0.0;
                if (logValues) {
                    viralR = betaK + sigma + Math.log(constantRTerm);
                } else {
                    viralR = betaK * sigma * constantRTerm;
                }
                sumRs += viralR;
                Rs.add(viralR);
                
            }
            
            // Compute mean and variance of betas
            double meanBeta = sumBetas / infections;
            double sumSqrDevsBeta = 0.0;
            for (int i = 0; i < infections; i++) {
                final double sqrDev = Math.pow((betas.get(i) - meanBeta),2);
                sumSqrDevsBeta += sqrDev;
            }
            final double varBeta = sumSqrDevsBeta / infections;
            fitDist.meanBeta = meanBeta;
            fitDist.varBeta = varBeta;

            // Compute mean and variance of sigmas
            double meanSigma = sumSigmas / infections;
            double sumSqrDevsSigmas = 0.0;
            for (int i = 0; i < infections; i++) {
                final double sqrDev = Math.pow((sigmas.get(i) - meanSigma),2);
                sumSqrDevsSigmas += sqrDev;
            }
            final double varSigma = sumSqrDevsSigmas / infections;
            fitDist.meanSigma = meanSigma;
            fitDist.varSigma = varSigma;
            
            // Compute mean and variance of Rs
            double meanR = sumRs / infections;
            double sumSqrDevsR = 0.0;
            for (int i = 0; i < infections; i++) {
                final double sqrDev = Math.pow((Rs.get(i) - meanR),2);
                sumSqrDevsR += sqrDev;
            }
            final double varR = sumSqrDevsR / infections;
            fitDist.meanR = meanR;
            fitDist.varR = varR;
            
            // Compute covariance of beta and sigma
            double sumProdDevs = 0.0;
            for (int i = 0; i < infections; i++) {
                final double prodDev = (betas.get(i) - meanBeta) * (sigmas.get(i) - meanSigma);
                sumProdDevs += prodDev;
            }
            final double covBetaSigma = sumProdDevs / infections;
            fitDist.covBetaSigma = covBetaSigma;
            
            return fitDist;
            
            
        }
        
        public ArrayList<Double> getMeanTypeImmunities(ArrayList<Integer> typeList) {
            
            ArrayList<Double> typeImmunities = new ArrayList<Double>(); //store average immunity in host population for each type            
            for (int t = 0; t < typeList.size(); t++) {
                
                int type = typeList.get(t);
                double meanTypeImmunity = 0;
                if (Parameters.hostImmuneHistorySampleCount >= getN()) {
                    
                    // If the number of host immune histories to sample is greater than the host pop size, loop through each host rather than sampling randomly
                    List<Host> allHosts = new ArrayList<Host>(); //includes S, I, and R hosts
                    allHosts.addAll(susceptibles); allHosts.addAll(infecteds); allHosts.addAll(recovereds);
                    double sumTypeImmunity = 0;
                    for (int s = 0; s < allHosts.size(); s++) {
                        Host sH = allHosts.get(s);
                        Phenotype[] history = sH.getHistory();
                        final double sigma = AntigenicTree.getClosestDistance(type, history); // Find closest phenotype in history in terms of logKb
                        sumTypeImmunity += sigma;  
                    }
                    meanTypeImmunity = sumTypeImmunity / Parameters.hostImmuneHistorySampleCount;
                    
                } else {
                    
                    // Sample (n = Parameters.hostImmmuneHistorySampleCount) host immune histories:
                    double sumTypeImmunity = 0;
                    for (int s = 0; s < Parameters.hostImmuneHistorySampleCount; s++) {
                        Host sH = getRandomHost();
                        Phenotype[] history = sH.getHistory();
                        final double sigma = AntigenicTree.getClosestDistance(type, history); // Find closest phenotype in history in terms of sigma
                        sumTypeImmunity += sigma; 
                    }
                    meanTypeImmunity = sumTypeImmunity / Parameters.hostImmuneHistorySampleCount;
                    
                }
                typeImmunities.add(meanTypeImmunity);
            
            }
            return typeImmunities;
                
        }
        
        // Not using this since absorbed into getViralFitnessDist
        public double getAntigenicLoadDistribution() {
            
            // New method for computing mean and variance of the antigenic load experienced by viral population
            // For each unique antigenic type, immunity is averaged over a set of sampled hosts' immune histories
            
            ArrayList<Integer> typeList = getAntigenicTypes(); // returns list of all unique antigenic types
            ArrayList<Integer> typeCounts = getAntigenicTypeCounts(typeList); // returns counts for each antigenic types ordered as in typeList
            ArrayList<Double> typeImmunities = new ArrayList<Double>(); //store average immunity in host population for each type
            double weightedSumOfTypes = 0; // the sum of the average immunities to each type weighted by typeCounts;
            
            for (int t = 0; t < typeList.size(); t++) {
                
                int type = typeList.get(t);
                double meanTypeImmunity = 0;
                
                if (Parameters.hostImmuneHistorySampleCount >= getN()) {
                    
                    // If the number of host immune histories to sample is greater than the host pop size, loop through each host rather than sampling randomly
                    List<Host> allHosts = new ArrayList<Host>(); //includes S, I, and R hosts
                    allHosts.addAll(susceptibles); allHosts.addAll(infecteds); allHosts.addAll(recovereds);
                    double sumTypeImmunity = 0;
                    for (int s = 0; s < allHosts.size(); s++) {
                        Host sH = allHosts.get(s);
                        Phenotype[] history = sH.getHistory();
                        final double sigma = AntigenicTree.getClosestDistance(type, history); // Find closest phenotype in history in terms of logKb
                        sumTypeImmunity += sigma;  
                    }
                    meanTypeImmunity = sumTypeImmunity / Parameters.hostImmuneHistorySampleCount;
                    
                } else {
                    
                    // Sample (n = Parameters.hostImmmuneHistorySampleCount) host immune histories:
                    double sumTypeImmunity = 0;
                    for (int s = 0; s < Parameters.hostImmuneHistorySampleCount; s++) {
                        Host sH = getRandomHost();
                        Phenotype[] history = sH.getHistory();
                        final double sigma = AntigenicTree.getClosestDistance(type, history); // Find closest phenotype in history in terms of sigma
                        sumTypeImmunity += sigma;  
                        // For Sean's antibody binding model
                        //final double logKbDelta = AntigenicTree.getClosestDistance(type, history); // Find closest phenotype in history in terms of logKb
                        //final double logKb = Parameters.logKbMax - logKbDelta;
                        //sumTypeImmunity += Math.exp(logKb);  
                    }
                    meanTypeImmunity = sumTypeImmunity / Parameters.hostImmuneHistorySampleCount;
                    
                }
                
                weightedSumOfTypes += meanTypeImmunity * typeCounts.get(t);
                typeImmunities.add(meanTypeImmunity);
                    
            }
            
            final double meanHostImmunity = weightedSumOfTypes / getI(); // average weighted by frequencies of each antigenic variant
            
            double weightedSumSqrDevs = 0.0;
            for (int t = 0; t < typeList.size(); t++) {
                final double sqrDev = Math.pow((typeImmunities.get(t) - meanHostImmunity),2);
                weightedSumSqrDevs += sqrDev * typeCounts.get(t);
            }
            final double varHostImmunity = weightedSumSqrDevs / getI(); 
            
            return meanHostImmunity;
            
        }
        
        
        public int getAntigenicCount() {
            ArrayList<Integer> typeList = new ArrayList<Integer>();
            if (getI() > 0) {
                for (int i = 0; i < getI(); i++) {
                    Host h = infecteds.get(i);
                    Virus v = h.getInfection();
                    Phenotype p = v.getPhenotype();
                    int type = p.antigenicType();
                    if (!typeList.contains(type)) {
                        typeList.add(type);
                    }
                }
            }
            int types = typeList.size();
            return types;
        }
        
        public ArrayList<Integer> getAntigenicTypes() {
            ArrayList<Integer> typeList = new ArrayList<Integer>();
            if (getI() > 0) {
                for (int i = 0; i < getI(); i++) {
                    Host h = infecteds.get(i);
                    Virus v = h.getInfection();
                    Phenotype p = v.getPhenotype();
                    int type = p.antigenicType();
                    if (!typeList.contains(type)) {
                        typeList.add(type);
                    }
                }
            }
            return typeList;
        }
        
        public ArrayList<Integer> getAntigenicTypeCounts(ArrayList<Integer> typeList) {
            ArrayList<Integer> typeCounts = new ArrayList<Integer>();
            for (int t = 0; t < typeList.size(); t++) {
                typeCounts.add(0);
            }
            if (getI() > 0) {
                for (int i = 0; i < getI(); i++) {
                    Host h = infecteds.get(i);
                    Virus v = h.getInfection();
                    Phenotype p = v.getPhenotype();
                    int type = p.antigenicType();
                    int typeIndex = typeList.indexOf(type);
                    typeCounts.set(typeIndex, typeCounts.get(typeIndex) + 1);
                    
                }
            }
            return typeCounts;
        }

        
        public ArrayList<Integer> getLoadDistribution() {
            
            ArrayList<Integer> dist = new ArrayList<Integer>();
            int classes = 50;
            int slot;
            for (int c = 0; c < classes; c++) {
                dist.add(0);
            }
            int infections = getI();
            if (infections > 0) {
                for (int i = 0; i < infections; i++) {
                    Host h = infecteds.get(i);
                    Virus v = h.getInfection();
                    Phenotype p = v.getPhenotype();
                    slot = p.mutLoad();
                    if (slot >= classes) {
                        slot = classes - 1;
                    }
                    dist.set(slot, dist.get(slot) + 1);
                }
            }
            
            return dist;
        }
	
	public void removeSusceptible(int i) {
		int lastIndex = getS() - 1;
		Host lastHost = susceptibles.get(lastIndex);
		susceptibles.set(i,lastHost);
		susceptibles.remove(lastIndex);
	}	
	public void removeInfected(int i) {
		int lastIndex = getI() - 1;
		Host lastHost = infecteds.get(lastIndex);
		infecteds.set(i,lastHost);
		infecteds.remove(lastIndex);
	}
	public void removeRecovered(int i) {
		int lastIndex = getR() - 1;
		Host lastHost = recovereds.get(lastIndex);
		recovereds.set(i,lastHost);
		recovereds.remove(lastIndex);
	}	
	
        // Trevor's method: mutations can occur any time during infection
	public void stepForward() {
	
	//	resetCases();
		if (Parameters.swapDemography) {
			swap();
		} else {
			grow();
			decline();
		}
		contact(); //<-mutations have to occur at transmission events
		recover();
		if (Parameters.transcendental) { 
			loseImmunity(); 
		}
		mutate();
		sample();
	
	}
        
        // My method: mutations only occur at transmission events
        public void stepForwardMutAtTrans() {
	
	//	resetCases();
		if (Parameters.swapDemography) {
			swap();
		} else {
			grow();
			decline();
		}
                infectionsNow = getI();
                newParents = new ArrayList<Host>();
                newInfections = new ArrayList<Host>();
                computeMeanLoad();
                getInfectionsByClass();
		//contactMutAtTrans(); //<-mutations have to occur at transmission events
                contactByMutClass();
                //recoverFix();
                recoverByClass();
		if (Parameters.transcendental) { 
			loseImmunity(); 
		}
                //if (Parameters.waningImmunity) {
                    //waneImmunity();
                //}
		sample();
                
                // Methods for cleaing up host immune histories (no long necessary)
                //if (Parameters.day > 0 & Parameters.day % 2000 == 0) {
                    //cleanUpHistory(); // clean up based on a threshold distance from all current antigenic types (fast, but does not preserve exact immune structure)
                    //cleanUpDetailedHistories(); // removes antigenic types from immune histories AND preserves exact immune structure (but slow)
                //}
	
	}
        
        public void getInfectionsByClass() {
            
            int infections = getI();
            int max = 0;
            ArrayList<Integer> kList = new ArrayList<Integer>();
            if (infections > 0) {
                for (int i = 0; i < infections; i++) {
                    Host h = infecteds.get(i);
                    Virus v = h.getInfection();
                    Phenotype p = v.getPhenotype();
                    int k = p.mutLoad();
                    if (k > max) {
                        max = k;
                    }
                    kList.add(k);
                }
            }
            
            infectionsByClass = new ArrayList<ArrayList<Host>>();
            for (int k = 0; k <= max; k++) {
                infectionsByClass.add(new ArrayList<Host>());
            }
            
            for (int i = 0; i < infections; i++) {
                infectionsByClass.get(kList.get(i)).add(infecteds.get(i));
            }
            
        }
        
	
	// draw a Poisson distributed number of births and add these hosts to the end of the population list
	public void grow() {
		double totalBirthRate = getN() * Parameters.birthRate;
		int births = Random.nextPoisson(totalBirthRate);
		for (int i = 0; i < births; i++) {
			Host h = new Host();
			susceptibles.add(h);
		}
	}
	
	// draw a Poisson distributed number of deaths and remove random hosts from the population list
	public void decline() {
		// deaths in susceptible class
		double totalDeathRate = getS() * Parameters.deathRate;
		int deaths = Random.nextPoisson(totalDeathRate);
		for (int i = 0; i < deaths; i++) {
			if (getS()>0) {
				int sndex = getRandomS();
				removeSusceptible(sndex);
			}
		}		
		// deaths in infectious class		
		totalDeathRate = getI() * Parameters.deathRate;
		deaths = Random.nextPoisson(totalDeathRate);
		for (int i = 0; i < deaths; i++) {
			if (getI()>0) {
				int index = getRandomI();
				removeInfected(index);
			}
		}		
	}
	
	// draw a Poisson distributed number of births and reset these individuals
	public void swap() {
		// draw random individuals from susceptible class
		double totalBirthRate = getS() * Parameters.birthRate;
		int births = Random.nextPoisson(totalBirthRate);
		for (int i = 0; i < births; i++) {
			if (getS()>0) {
				int index = getRandomS();
				Host h = susceptibles.get(index);
				h.reset();
			}
		}		
		// draw random individuals from infected class
		totalBirthRate = getI() * Parameters.birthRate;
		births = Random.nextPoisson(totalBirthRate);
		for (int i = 0; i < births; i++) {
			if (getI()>0) {
				int index = getRandomI();
				Host h = infecteds.get(index);
				h.reset();
				removeInfected(index);
				susceptibles.add(h);
			}
		}	
		// draw random individuals from recovered class
		totalBirthRate = getR() * Parameters.birthRate;
		births = Random.nextPoisson(totalBirthRate);
		for (int i = 0; i < births; i++) {
			if (getR()>0) {
				int index = getRandomR();
				Host h = recovereds.get(index);
				h.reset();
				removeRecovered(index);
				susceptibles.add(h);
			}
		}			
	}

	// draw a Poisson distributed number of contacts and move from S->I based upon this
	public void contact() {

		// each infected makes I->S contacts on a per-day rate of beta * S/N
		double totalContactRate = getI() * getPrS() * Parameters.beta * Parameters.getSeasonality(deme);
		int contacts = Random.nextPoisson(totalContactRate);		
		for (int i = 0; i < contacts; i++) {
			if (getS()>0 && getI()>0) {
		
				// get indices and objects
				int index = getRandomI();
				int sndex = getRandomS();			
				Host iH = infecteds.get(index);			
				Host sH = susceptibles.get(sndex);						
				Virus v = iH.getInfection();
								
				// attempt infection
				Phenotype p = v.getPhenotype();		
				Phenotype[] history = sH.getHistory();
				double chanceOfSuccess = p.riskOfInfection(history);
				if (Random.nextBoolean(chanceOfSuccess)) {
					sH.infect(v,deme);
					removeSusceptible(sndex);
					infecteds.add(sH);
					cases++;
				}
			
			}
		}		
		
	}
        
        // draw a Poisson distributed number of contacts and move from S->I based upon this
	public void contactMutAtTrans() {

		// each infected makes I->S contacts on a per-day rate of beta * S/N
		double totalContactRate = (getI() + Parameters.externalMigration) * getPrS() * Parameters.beta * Parameters.getSeasonality(deme);
		int contacts = Random.nextPoisson(totalContactRate);
                //double varScaler = 0.05;
                //int contacts = (int) Math.round(Random.nextNormal(totalContactRate, Math.sqrt(varScaler * totalContactRate)));
		for (int i = 0; i < contacts; i++) {
			if (getS()>0 && getI()>0) {
		
				// get indices and objects
				int index = getRandomI();
                                while (newParents.contains(infecteds.get(index)) || newInfections.contains(infecteds.get(index))) { //or new infections?
                                    //System.out.println("Already a parent");
                                    index = getRandomI();
                                }
				int sndex = getRandomS();			
				Host iH = infecteds.get(index);			
				Host sH = susceptibles.get(sndex);						
				Virus v = iH.getInfection();
                                
				// attempt infection
				Phenotype p = v.getPhenotype();		
				Phenotype[] history = sH.getHistory();
				double chanceOfSuccess = p.riskOfInfection(history);
				if (Random.nextBoolean(chanceOfSuccess)) {
                                        // add mutations before transmission
                                        Virus mutV = v.mutate(); //but will this mutate parents virus? should leave intact
                                        Phenotype mutP = mutV.getPhenotype();
                                        if (!mutP.lethal()) {
                                            //sH.infect(mutV,deme);
                                            iH.infect(v, deme);
                                            sH.infectWithMutant(mutV, deme);
                                            removeSusceptible(sndex);
                                            infecteds.add(sH);
                                            cases++;
                                            newParents.add(iH);
                                            newInfections.add(sH);
                                        }
				}
			
			}
		}		
		
	}
        
        // draw a Poisson distributed number of contacts and move from S->I based upon this
	public void contactByMutClass() {
                
                int maxK = infectionsByClass.size();
                for (int k = 0; k < maxK; k++) {
                    
                    int infectionsK = infectionsByClass.get(k).size();
                    double betaK = Parameters.beta * Math.pow((1-Parameters.mutCost), k);
                    double propK = infectionsK / getI();

                    double totalContactRate = (infectionsK + (Parameters.externalMigration * propK)) * getPrS() * betaK * Parameters.getSeasonality(deme);
                    int contacts = (int) Math.round(Random.nextNormal(totalContactRate, Math.sqrt(Parameters.demoNoiseScaler * totalContactRate)));
                    
                    if (contacts > infectionsK) {
                        System.out.println("More contacts than infections in class!");
                        contacts = infectionsK;
                    }
                    
                    for (int i = 0; i < contacts; i++) {
                            if (getS()>0 && getI()>0) {

                                    // get indices and objects
                                    int randomInK = (int) Math.floor(Math.random() * infectionsK);
                                    Host testH = infectionsByClass.get(k).get(randomInK);
                                    while (newParents.contains(testH) || newInfections.contains(testH)) { //or new infections?
                                        //index = getRandomI();
                                        randomInK = (int) Math.floor(Math.random() * infectionsK);
                                        testH = infectionsByClass.get(k).get(randomInK);
                                    }
                                    int sndex = getRandomS();			
                                    Host iH = testH; //infecteds.get(index);			
                                    Host sH = susceptibles.get(sndex);						
                                    Virus v = iH.getInfection();

                                    // attempt infection
                                    Phenotype p = v.getPhenotype();		
                                    Phenotype[] history = sH.getHistory();
                                    double chanceOfSuccess = p.riskOfInfection(history);
                                    if (Random.nextBoolean(chanceOfSuccess)) {
                                            // add mutations before transmission
                                            Virus mutV = v.mutate();
                                            Phenotype mutP = mutV.getPhenotype();
                                            if (!mutP.lethal()) {
                                                //sH.infect(mutV,deme);
                                                iH.infect(v, deme);
                                                sH.infectWithMutant(mutV, deme);
                                                removeSusceptible(sndex);
                                                infecteds.add(sH);
                                                cases++;
                                                newParents.add(iH);
                                                newInfections.add(sH);
                                            }
                                    }

                            }
                    }
                }
		
	}
	
	// draw a Poisson distributed number of contacts and move from S->I based upon this
	// this deme is susceptibles and other deme is infecteds
	public void betweenDemeContact(HostPopulation hp) {

		// each infected makes I->S contacts on a per-day rate of beta * S/N
		double totalContactRate = hp.getI() * getPrS() * Parameters.beta * Parameters.betweenDemePro * Parameters.getSeasonality(deme);
		int contacts = Random.nextPoisson(totalContactRate);
		for (int i = 0; i < contacts; i++) {
			if (getS()>0 && hp.getI()>0) {
		
				// get indices and objects
				Host iH = hp.getRandomHostI();
				int sndex = getRandomS();
				Host sH = susceptibles.get(sndex);			
				Virus v = iH.getInfection();
				
				// attempt infection
				Phenotype p = v.getPhenotype();
				Phenotype[] history = sH.getHistory();
				double chanceOfSuccess = p.riskOfInfection(history);
				if (Random.nextBoolean(chanceOfSuccess)) {
					sH.infect(v,deme);
					removeSusceptible(sndex);
					infecteds.add(sH);
					cases++;
				}
			
			}
		}		
		
	}	
	
	// draw a Poisson distributed number of recoveries and move from I->S based upon this
	public void recover() {
		// each infected recovers at a per-day rate of nu
		double totalRecoveryRate = getI() * Parameters.nu;
		int recoveries = Random.nextPoisson(totalRecoveryRate);
		for (int i = 0; i < recoveries; i++) {
			if (getI()>0) {
				int index = getRandomI();
				Host h = infecteds.get(index);
				h.clearInfection();
				removeInfected(index);
				if (Parameters.transcendental) {
					recovereds.add(h);
				} else {
					susceptibles.add(h);
				}

			}
		}			
	}
        
        // draw a Poisson distributed number of recoveries and move from I->S based upon this
	public void recoverFix() {
		// each infected recovers at a per-day rate of nu
		double totalRecoveryRate = infectionsNow * Parameters.nu;
		int recoveries = Random.nextPoisson(totalRecoveryRate);
                //double varScaler = 0.05;
                //int recoveries = (int) Math.round(Random.nextNormal(totalRecoveryRate, Math.sqrt(varScaler * totalRecoveryRate)));
		for (int i = 0; i < recoveries; i++) {
			if (infectionsNow>0) {
                                
				int index = getRandomI();
                                //while (newInfections.contains(infecteds.get(index))) {
                                  //  index = getRandomI();
                                //}
				Host h = infecteds.get(index);
				h.clearInfection();
				removeInfected(index);
				if (Parameters.transcendental) {
					recovereds.add(h);
				} else {
					susceptibles.add(h);
				}

			}
		}			
	}
        
        // draw a Poisson distributed number of recoveries and move from I->S based upon this
	public void recoverByClass() {
            
                int maxK = infectionsByClass.size();
                for (int k = 0; k < maxK; k++) {
                    
                    int infectionsK = infectionsByClass.get(k).size();
                
                    // each infected recovers at a per-day rate of nu
                    double totalRecoveryRate = infectionsK * Parameters.nu;
                    
                    int recoveries = (int) Math.round(Random.nextNormal(totalRecoveryRate, Math.sqrt(Parameters.demoNoiseScaler * totalRecoveryRate)));
                    
                    if (recoveries > infectionsK) {
                        System.out.println("More recoveries than infections in class!");
                        recoveries = infectionsK;
                    }
                    
                    for (int i = 0; i < recoveries; i++) {
			if (infectionsNow>0) {
                                
                                int randomInK = (int) Math.floor(Math.random() * infectionsByClass.get(k).size());
                                Host h = infectionsByClass.get(k).get(randomInK);
                                int indexInInfecteds = infecteds.indexOf(h);
				h.clearInfection();
                                infectionsByClass.get(k).remove(randomInK);
				removeInfected(indexInInfecteds);
				if (Parameters.transcendental) {
					recovereds.add(h);
				} else {
					susceptibles.add(h);
				}

			}
                    }
                }
	}

	
	// draw a Poisson distributed number of R->S 
	public void loseImmunity() {
		// each recovered regains immunity at a per-day rate
		double totalReturnRate = getR() * Parameters.immunityLoss;
		int returns = Random.nextPoisson(totalReturnRate);
		for (int i = 0; i < returns; i++) {
			if (getR()>0) {
				int index = getRandomR();
				Host h = recovereds.get(index);
                                h.reset(); //Added this for waning immunity
				removeRecovered(index);
				susceptibles.add(h);
			}
		}			
	}	
	
	// draw a Poisson distributed number of mutations and mutate based upon this
	// mutate should not impact other Virus's Phenotypes through reference
	public void mutate() {
		// each infected mutates at a per-day rate of mu
		double totalMutationRate = getI() * Parameters.muPhenotype;
		int mutations = Random.nextPoisson(totalMutationRate);
		for (int i = 0; i < mutations; i++) {
			if (getI()>0) {
				int index = getRandomI();
				Host h = infecteds.get(index);
				h.mutate();
			}
		}			
	}	
	
	// draw a Poisson distributed number of samples and add them to the VirusSample
	// only sample after burnin is completed
	public void sample() {
		if (getI()>0 && Parameters.day >= Parameters.burnin) {
		
			double totalSamplingRate = Parameters.tipSamplingRate;
			if (Parameters.tipSamplingProportional) {
				totalSamplingRate *= getI();
			} 
			
			int samples = Random.nextPoisson(totalSamplingRate);
			
			if(samples < 10){
				samples = 10;
			}
			if(getI() < 10 || samples > getI() ){
				samples = getI();
			}
			for (int i = 0; i < samples; i++) {
				int index = getRandomI();
				Host h = infecteds.get(index);
				
				if(!h.getHasBeenSampled()){
					Virus v = h.getInfection();
					VirusTree.add(v);
					h.markSampled();
				}
				
			}	
		}
	}
		
	// through current infected population assigning ancestry as trunk
	public void makeTrunk() {
		for (int i = 0; i < getI(); i++) {
                        //int i = getRandomI();
			Host h = infecteds.get(i);
			Virus v = h.getInfection();
			v.makeTrunk();
			while (v.getParent() != null) {
				v = v.getParent();
				if (v.isTrunk()) {
					break;
				} else {
					v.makeTrunk();
				}
			}
		}
	}	
	
	public void updateDiversity() {

		diversity = 0.0;
		tmrca = 0.0;
		antigenicDiversity = 0.0;		
		netau = 0.0;
		serialInterval = 0.0;
		
		if (getI()>1) { 
		
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
			antigenicDiversity /= (double) sampleCount;		
			netau = coalOpp / coalCount;
			serialInterval /= (double) sampleCount;
		
		}
		
	}
        
        public void printMutations(PrintStream stream) {
		ArrayList<Integer> dist = getLoadDistribution();
                for (int i = 0; i < dist.size(); i++) {
                    //String out = Integer.toString(dist.get(i));
                    stream.printf("%d", dist.get(i));
                    if (i < dist.size() - 1) {
                        stream.printf("%s", ", ");
                    }
                }
		//stream.printf("\t%.4f\t%.4f\t%.4f\t%.5f\t%.4f\t%d\t%d\t%d\t%d\t%d", getDiversity(), getTmrca(), getNetau(), getSerialInterval(), getAntigenicDiversity(), getN(), getS(), getI(), getR(), getCases());
	}
		
	public void printState(PrintStream stream) {
		updateDiversity();
		stream.printf("\t%.4f\t%.4f\t%.4f\t%.5f\t%.4f\t%d\t%d\t%d\t%d\t%d", getDiversity(), getTmrca(), getNetau(), getSerialInterval(), getAntigenicDiversity(), getN(), getS(), getI(), getR(), getCases());
	}	
	
	public void printHeader(PrintStream stream) {
		stream.printf("\t%sDiversity\t%sTmrca\t%sNetau\t%sSerialInterval\t%sAntigenicDiversity\t%sN\t%sS\t%sI\t%sR\t%sCases", name, name, name, name, name, name, name, name, name, name);
	}
	
	// reset population to factory condition
	public void reset() {
	
		// clearing lists
		susceptibles.clear();
		infecteds.clear();
		recovereds.clear();
		
		int initialR = 0;
		if (Parameters.transcendental) {
			initialR = (int) ((double) Parameters.initialNs[deme] * Parameters.initialPrT);
		}
	
		// fill population with susceptibles
		int initialS = Parameters.initialNs[deme] - Parameters.initialI - initialR;
		for (int i = 0; i < initialS; i++) {
			Host h = new Host();			
			susceptibles.add(h);
		}
		
		// fill population with recovereds
		for (int i = 0; i < initialR; i++) {
			Host h = new Host();			
			recovereds.add(h);
		}		
		
		if (deme == 1) {
		
			// infect some individuals
			for (int i = 0; i < 3*Parameters.initialI; i++) {
				Virus v = new Virus(Parameters.urVirus, deme);
				Host h = new Host(v);
				infecteds.add(h);
			}	
		
		}
		
	}
	
	public void printHostPopulation(PrintStream stream) {
		
		// step through susceptibles and print
		for (int i = 0; i < getS(); i++) {
			Host h = susceptibles.get(i);
			stream.print(deme + ":");
			h.printInfection(stream);
			stream.print(":");
			h.printHistory(stream);
			stream.println();
		}
		
		// step through infecteds and print
		for (int i = 0; i < getI(); i++) {
			Host h = infecteds.get(i);
			stream.print(deme + ":");
			h.printInfection(stream);
			stream.print(":");
			h.printHistory(stream);
			stream.println();
		}
		
		// step through recovereds and print
		for (int i = 0; i < getR(); i++) {
			Host h = recovereds.get(i);
			stream.print(deme + ":");
			h.printInfection(stream);
			stream.print(":");
			h.printHistory(stream);
			stream.println();
		}		
	
	}
        
        /**
         * Old methods that should no longer be used
         */
        
        // Defunct method!!!
        public void stepForwardMutAtTransConstantSize() {
	
	//	resetCases();
                int infectionsBeforeDeaths = getI();
		if (Parameters.swapDemography) {
			swap();
		} else {
			grow();
			decline();
		}
                infectionsNow = getI();
                newParents = new ArrayList<Host>();
                newInfections = new ArrayList<Host>();
                computeMeanLoad();
		//contactMutAtTrans(); //<-mutations have to occur at transmission events
                getInfectionsByClass();
                contactByMutClass();
                int infectionsAfterBirths = getI();
                recoveriesNow = infectionsAfterBirths - infectionsBeforeDeaths;
                recoverByClassConstantSize();
		if (Parameters.transcendental) { 
			loseImmunity(); 
		}
                //if (Parameters.day >= Parameters.tipSamplingStartDay) {
                    sample();
                //}
	
	}
        
        // Defunct
        public void recoverByClassConstantSize() {
            
                int maxK = infectionsByClass.size();
                double infNowDouble = (double) infectionsNow;
                
                ArrayList<Integer> recoveriesK = new ArrayList<Integer>();
                int totalR = 0;
                for (int k = 0; k < maxK; k++) {
                    int infectionsK = infectionsByClass.get(k).size();
                    double infKDouble = (double) infectionsK;
                    int recoveries = (int) Math.round((infKDouble / infNowDouble) * recoveriesNow);
                    recoveriesK.add(recoveries);
                    totalR += recoveries;
                }
                int remainder = recoveriesNow - totalR;
                if (remainder < 0) { //too many
                    remainder = Math.abs(remainder);
                    for (int i = 0; i < remainder; i++) {
                        int randLoc = (int) Math.floor(Math.random() * maxK);
                        while (recoveriesK.get(randLoc) <= 1) {
                            randLoc = (int) Math.floor(Math.random() * maxK);
                        }
                        recoveriesK.set(randLoc, recoveriesK.get(randLoc) - 1);
                        totalR--;
                    }
                }
                if (remainder > 0) { //too few
                    for (int i = 0; i < remainder; i++) {
                        int randLoc = (int) Math.floor(Math.random() * maxK);
                        while (recoveriesK.get(randLoc) >= infectionsByClass.get(randLoc).size()) {
                            randLoc = (int) Math.floor(Math.random() * maxK);
                        }
                        recoveriesK.set(randLoc, recoveriesK.get(randLoc) + 1);
                        totalR++;
                    }
                }
                //System.out.println();
                
                
                for (int k = 0; k < maxK; k++) {
                    
                    int infectionsK = infectionsByClass.get(k).size();
                
                    // each infected recovers at a per-day rate of nu
                    //double totalRecoveryRate = infectionsK * Parameters.nu;
                    
                    //double varScaler = 0.05;
                    //double oldSD = Math.sqrt(totalRecoveryRate);
                    //double newSD = Math.sqrt(varScaler * totalRecoveryRate);
                    //if (newSD > oldSD) {
                        //System.out.println("Increase in removal variance");
                    //}
                    //int recoveries = (int) Math.round(Random.nextNormal(totalRecoveryRate, Math.sqrt(Parameters.demoNoiseScaler * totalRecoveryRate)));
                    
                    //double infKDouble = (double) infectionsK;
                    //int recoveries = (int) Math.round((infKDouble / infNowDouble) * recoveriesNow);
                    
                    int recoveries = recoveriesK.get(k);
                    
                    if (recoveries > infectionsK) {
                        System.out.println("More recoveries than infections in class!");
                        recoveries = infectionsK;
                    }
                    
                    for (int i = 0; i < recoveries; i++) {
			if (infectionsNow>0) {
                                
                                int randomInK = (int) Math.floor(Math.random() * infectionsByClass.get(k).size());
                                Host h = infectionsByClass.get(k).get(randomInK);
                                int indexInInfecteds = infecteds.indexOf(h);
				h.clearInfection();
                                infectionsByClass.get(k).remove(randomInK);
				removeInfected(indexInInfecteds);
				if (Parameters.transcendental) {
					recovereds.add(h);
				} else {
					susceptibles.add(h);
				}

			}
                    }
                }
	}
        
        //Defunct
        // draw a Poisson distributed number of recoveries and move from I->S based upon this
	public void recoverConstantSize() {
		// each infected recovers at a per-day rate of nu
		//double totalRecoveryRate = infectionsNow * Parameters.nu;
		int recoveries = recoveriesNow; //Random.nextPoisson(totalRecoveryRate);
		for (int i = 0; i < recoveries; i++) {
			if (infectionsNow>0) {
                                
				int index = getRandomI();
                                //while (newInfections.contains(infecteds.get(index))) {
                                    //index = getRandomI();
                                //}
				Host h = infecteds.get(index);
				h.clearInfection();
				removeInfected(index);
				if (Parameters.transcendental) {
					recovereds.add(h);
				} else {
					susceptibles.add(h);
				}

			}
		}			
	}
        
        // Defunct
        public void cleanUpDetailedHistories() {
            
            ArrayList<Integer> presentTypes = getAntigenicTypes();
            HashSet<Integer> retainTypes = new HashSet<Integer>();
            double distance; double min;
            
            for (int sndex = 0; sndex < getS(); sndex++) {
                Host sH = susceptibles.get(sndex);
                ArrayList<Integer> hArray = sH.getHistoryArray();
                if (hArray.size() > 1) {
                    HashSet<Integer> retain = new HashSet<Integer>();
                    for (int t : presentTypes) {
                        min = Double.POSITIVE_INFINITY;
                        Integer hold = 0;
                        for (int s : hArray) {
                            distance = AntigenicTree.getDistance(s, t);
                            if (distance < min) {
                                min = distance;
                                hold = s;
                            }
                        }
                        if (min < 1.0) {    // added this for no sibling immunity case
                            retain.add(hold);
                        }
                    }
                    retainTypes.addAll(retain);
                    hArray.removeAll(retain);
                    if (hArray.size() > 0) {
                        sH.removeFromHistory(hArray);
                    }
                } else if (hArray.size() == 1) {
                    retainTypes.addAll(hArray);
                }   
            }
            
            for (int index = 0; index < getI(); index++) {
                Host iH = infecteds.get(index);
                ArrayList<Integer> hArray = iH.getHistoryArray();
                if (hArray.size() > 1) {
                    HashSet<Integer> retain = new HashSet<Integer>();
                    for (int t : presentTypes) {
                        min = Double.POSITIVE_INFINITY;
                        Integer hold = 0;
                        for (int s : hArray) {
                            distance = AntigenicTree.getDistance(s, t);
                            if (distance < min) {
                                min = distance;
                                hold = s;
                            }
                        }
                        if (min < 1.0) {    // addes this for no sibling immunity case
                            retain.add(hold);
                        }
                    }
                    retainTypes.addAll(retain);
                    hArray.removeAll(retain);
                    if (hArray.size() > 0) {
                        iH.removeFromHistory(hArray);
                    }
                } else if (hArray.size() == 1) {
                    retainTypes.addAll(hArray);
                }
            }
            
            ArrayList<Integer> removeTypes = new ArrayList<Integer>();
            removeTypes.addAll(AntigenicTree.getTreeTypes());
            removeTypes.removeAll(retainTypes);
            removeTypes.removeAll(presentTypes); //Not sure if need this
            for (int type : removeTypes) {
                AntigenicTree.remove(type);
            }
            System.out.println("Removed " + removeTypes.size() + " from population's antigenic history");
            
        }
        
        
        // Defunct
        public void cleanUpHistory() {
            
            ArrayList<Integer> removeList = new ArrayList<Integer>();
            ArrayList<Integer> presentTypes = getAntigenicTypes();
            int numPresent = presentTypes.size();
            ArrayList<Integer> treeTypes = AntigenicTree.getTreeTypes(); //return currTypes
            int removed = 0;
            for (int i = 0; i < treeTypes.size(); i++) {
                boolean retain = false;
                int counter = 0;
                int type;
                double distance;
                while (!retain & counter < numPresent) {
                    type = presentTypes.get(counter);
                    distance = AntigenicTree.getDistance(treeTypes.get(i),type);
                    //System.out.println("Distance = " + distance);
                    if (distance <= Parameters.cleanUpDistance) {
                        retain = true;
                    }
                    counter++;
                }
                if (!retain) {
                    //System.out.println(treeTypes.get(i));
                    AntigenicTree.remove(treeTypes.get(i));
                    removeList.add(treeTypes.get(i));
                    removed++;
                }
            }
            
            if (removeList.size() > 0) {
                for (int sndex = 0; sndex < getS(); sndex++) {
                    Host sH = susceptibles.get(sndex);
                    sH.removeFromHistory(removeList);
                }
                for (int index = 0; index < getI(); index++) {
                    Host iH = infecteds.get(index);
                    iH.removeFromHistory(removeList);
                }
                //Have to clean up history for recoverds too if transcendental immunity
            }
            
            System.out.println("Removed " + removed + " from population's antigenic history");
            //System.out.println(removeList);
        }
				
}
