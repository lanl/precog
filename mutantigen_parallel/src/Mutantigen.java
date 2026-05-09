/* Implements an individual-based model in which the infection's genealogical history is tracked through time */

class Mutantigen {

    public static String globalSimNum;
    
    public static void main(String[] args) {
    
    		globalSimNum = args[1];

		// initialize random number generator
		cern.jet.random.AbstractDistribution.makeDefaultGenerator();
		
		// initialize static parameters
		Parameters.load(args[0]);		
		Parameters.initialize();
                
                // initialize antigenic tree
                AntigenicTree.initialize();
		
		// run simulation
		Simulation sim = new Simulation(args[1]);
		sim.run();	
		
	}
   	
}
