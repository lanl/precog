/* A human individual that harbors viruses and immunity */

import java.util.*;
import java.io.*;
import java.util.regex.*;

public class Host {

	// fields
	private Virus infection;												
	private Phenotype[] immuneHistory = new Phenotype[0];
	private boolean hasBeenSampled;
	
	// naive host
	public Host() {
		hasBeenSampled = false;
		initializeHistory();		
	}
	
	// initial infected host
	public Host(Virus v) {
		hasBeenSampled = false;
		infection = v;
		initializeHistory();
	}
	
	// checkpointed host
	public Host(int d, String sVirus, String sHist) {
		if (!sVirus.equals("n")) {
			Pattern rc = Pattern.compile(",");
    		String[] traitList = rc.split(sVirus);
    		double x = Double.parseDouble(traitList[0]);
    		double y = Double.parseDouble(traitList[1]);
			Phenotype p = PhenotypeFactory.makeArbitaryPhenotype(x,y);
			infection = new Virus(Parameters.urVirus,d,p);
		}
		if (!sHist.equals("n")) {
			Pattern rsc = Pattern.compile(";");
    		String[] phenotypeList = rsc.split(sHist);
    		for (int i = 0; i < phenotypeList.length; i++) {
				Pattern rc = Pattern.compile(",");
    			String[] traitList = rc.split(phenotypeList[i]);
    			double x = Double.parseDouble(traitList[0]);
    			double y = Double.parseDouble(traitList[1]);
				Phenotype p = PhenotypeFactory.makeArbitaryPhenotype(x,y);
				addToHistory(p);
			}
		}		
	}

	public void markSampled(){
		hasBeenSampled = true;
	}

	// sometimes start with immunity	
	public void initializeHistory() {
		double chanceOfSuccess = Parameters.initialPrR;
                if (Parameters.backgroundImmunity) {
                    if (Random.nextBoolean(chanceOfSuccess)) {	
			Phenotype p = Parameters.urImmunity;    //now antigenic strain 1, not 0
			addToHistory(p);
                    } else {
                        Phenotype p = PhenotypeFactory.makeBackgroundPhenotype();
                        addToHistory(p);
                    }
                    
                } else {
                    if (Random.nextBoolean(chanceOfSuccess)) {	
			Phenotype p = Parameters.urImmunity;
			addToHistory(p);
                    }
		}	
	}
	
	public void addToHistory(Phenotype p) {
		Phenotype[] newHistory = new Phenotype[immuneHistory.length + 1];
		for (int i = 0; i < immuneHistory.length; i++) {
			newHistory[i] = immuneHistory[i];
		}
		newHistory[immuneHistory.length] = p;
		immuneHistory = newHistory;
	}
	
	// infection methods
	public void reset() {
		infection = null;
		immuneHistory = new Phenotype[0];
	}
        
        public void removeFromHistory(ArrayList<Integer> removeList) {
            ArrayList<Integer> retain = new ArrayList<Integer>();
            for (int i = 0; i < immuneHistory.length; i++) {
                if (!removeList.contains(immuneHistory[i].antigenicType())) {
                    retain.add(i);
                }
            }
            Phenotype[] newHistory = new Phenotype[retain.size()];
            for (int i = 0; i < retain.size(); i++) {
                newHistory[i] = immuneHistory[retain.get(i)];
            }
            immuneHistory = newHistory;
        }
	
	public boolean isInfected() {
		boolean infected = false;
		if (infection != null) {
			infected = true;
		}
		return infected;
	}

	public boolean getHasBeenSampled() {
		return hasBeenSampled;
	}

	public Virus getInfection() {
		return infection;
	}
	public void infect(Virus pV, int d) {
		Virus nV = new Virus(pV, d);
		infection = nV;
	}
        public void infectWithMutant(Virus pV, int d) {
		//Virus nV = new Virus(pV, d);
		infection = pV;
	}
	public void clearInfection() {
		Phenotype p = infection.getPhenotype();
		addToHistory(p);
		infection = null;
		hasBeenSampled = false;
	}
	public int getHistoryLength() {
		return immuneHistory.length;
	}
	
	// make a new virus with the mutated phenotype
	public void mutate() {
		Virus mutV = infection.mutate();
		infection = mutV;
	}
	
	// history methods
	public Phenotype[] getHistory() {
		return immuneHistory;
	}
        
        public ArrayList<Integer> getHistoryArray() {
                ArrayList<Integer> iList = new ArrayList<Integer>();
            	for (int i = 0; i < immuneHistory.length; i++) {
                        iList.add(immuneHistory[i].antigenicType());
		}
                return iList;
        }
	
	public void printHistory() {
		for (int i = 0; i < immuneHistory.length; i++) {
			System.out.println(immuneHistory[i]);
		}
	}

	public void printInfection(PrintStream stream) {
		if (infection != null) {
			stream.print(infection.getPhenotype());
		}
		else {
			stream.print("n");
		}
	}
	
	public void printHistory(PrintStream stream) {
		if (immuneHistory.length > 0) {
			stream.print(immuneHistory[0]);
			for (int i = 1; i < immuneHistory.length; i++) {
				stream.print(";" + immuneHistory[i]);
			}
		}
		else {
			stream.print("n");
		}
	}	
		
	public String toString() {
		return Integer.toHexString(this.hashCode());
	}	
	
}
