/* Antigenic phenotype present in individual Viruses and within Hosts as immune history */
/* Should be able to calculate distance and cross-immunity between two phenotypes */
/* Moving up to multiple dimensions is non-trivial and requires thought on the implementation */
/* Multiple Viruses can reference a single Phenotype object */

import static java.lang.Math.*;
import java.util.*;

public class MutLoadPhenotype implements Phenotype {

	// fields
    
        //Don't need these
	private double traitA;
	private double traitB;
        
        //Mutation parameters
        private int mutLoad;
        private boolean lethal;
        
        //Antigenic type
        private int antigenType;
	
	// constructor
	public MutLoadPhenotype() {
	
	}
        
        public MutLoadPhenotype(int load, boolean lethality, int type) {
            mutLoad = load;
            lethal = lethality;
            antigenType = type;
        }
        
        public MutLoadPhenotype(double theta) {
            mutLoad = Random.nextPoisson(theta);
            lethal = false;
            antigenType = 0;
        }
        
        public double riskOfInfection( Phenotype[] history) {
	
		// find closest phenotype in history
		double closestDistance = 100.0;
		if (history.length > 0) {
			for (int i = 0; i < history.length; i++) {
				double thisDistance = distance(history[i]);
				if (thisDistance < closestDistance) {
					closestDistance = thisDistance;
				}
			}
		} 
		
		double antigenicRisk =  closestDistance;
                if (antigenicRisk > 1.0) {
                    antigenicRisk = 1.0;
                }
                
                double risk = 1.0;
                if (lethal) {
                    risk = 0.0;
                } //else {
                    //risk = (Parameters.beta * Math.pow((1-Parameters.mutCost),mutLoad)) / Parameters.beta; //standardized risk
                //}
                //System.out.println();
                
                risk *= antigenicRisk;
				
		return risk;
		
	}
       
       	// raw antigenic distance between two phenotypes
	public double distance(Phenotype p) {
                MutLoadPhenotype p2 = (MutLoadPhenotype) p;
                /**
                 * Distances could be precomputed between all types!!
                 */
                double dist = AntigenicTree.getDistance(antigenType, p2.antigenType);
                //System.out.println("Distance = " + dist);
                return dist;
	}
       
       	public Phenotype mutate() {
		
		// direction of mutation
		//double theta = 0;
		//if (Parameters.mut2D) {
		//	theta = Random.nextDouble(0,2*Math.PI);
		//} else {
		//	if (Random.nextBoolean(0.5)) { theta = 0; }
		//	else { theta = Math.PI; }
		//}
		
		// size of mutation
		//double alpha = (Parameters.meanStep *  Parameters.meanStep) / (Parameters.sdStep * Parameters.sdStep);
		//double beta = (Parameters.sdStep * Parameters.sdStep) / Parameters.meanStep;
		//double r = Random.nextGamma(alpha, beta);
			
		//double mutA = getTraitA() + r * Math.cos(theta);
		//double mutB = getTraitB() + r * Math.sin(theta);
		//Phenotype mutP = new GeometricPhenotype(mutA,mutB);
                
                //int mutations = Random.nextPoisson(Parameters.lambda);
                //int newLoad = mutLoad + mutations;
                //double probNoLethals = Math.pow((1-Parameters.probLethal), mutations);
                //boolean lethality = false;
                //if (probNoLethals < Math.random()) {
                    //lethality = true;
                //}
                
                int mutations = Random.nextPoisson(Parameters.lambda);
                int newLoad = mutLoad;
                boolean lethality = false;
                //double prBenificial = Parameters.epsilonSlope * Parameters.meanLoad;
                for (int m = 0; m < mutations; m++) {
                    if (Math.random() < Parameters.probLethal) {
                        lethality = true;
                    } else {
                        if (Math.random() < Parameters.epsilon) {
                            if (newLoad > 0) {newLoad--;}
                        } else {
                            newLoad++;    
                        }
                    }
                }
                
                int type = antigenType;
                if (Parameters.day >= Parameters.antigenicEvoStartDay) {
                    int antigenicMutations = Random.nextPoisson(Parameters.lambdaAntigenic);
                    if (antigenicMutations > 0 & !lethality) {
                        double distanceFromParent = 0.0;
                        for (int m = 0; m < antigenicMutations; m++) {
                            if (Parameters.antigenicGammaShape == 1.0) {
                                distanceFromParent += Random.nextExponential(Parameters.meanAntigenicSize);
                            } else {
                                distanceFromParent += Random.nextGamma(Parameters.antigenicGammaShape, (Parameters.meanAntigenicSize / Parameters.antigenicGammaShape));
                            }
                        }
                        //System.out.println("New ant mutant distance = " + distanceFromParent);
                        if (distanceFromParent >= Parameters.thresholdAntigenicSize) {
                            type = AntigenicTree.add(antigenType, distanceFromParent);
                        }
                    }
                }
                
                Phenotype mutP = new MutLoadPhenotype(newLoad, lethality, type);
                
		return mutP;
				
	}
        
        public String toString() {
		String fullString = String.format("%.4f,%.4f", traitA, traitB);
		return fullString;
	}
        
        public int mutLoad() {
            return mutLoad;
        }
        
        public boolean lethal() {
            return lethal;
        }
        
        public int antigenicType() {
            return antigenType;
        }
       
       
}
