/* Stores a list of Viruses that have sampled during the course of the simulation */

import java.util.*;
import java.io.*;
import static java.lang.Math.*;
import cern.colt.list.DoubleArrayList;

public class VirusTree {

	// fields
	private static Virus root = Parameters.urVirus;	
	private static List<Virus> tips = new ArrayList<Virus>();
        
        private static ArrayList<Virus> samples;
	
	public static double xMin;
	public static double xMax;
	public static double yMin;
	public static double yMax;
	public static double zMin;
	public static double zMax;
        
        public static ArrayList<Integer> tipTypes;
        public static ArrayList<Integer> majorTypes;
        public static ArrayList<Integer> firstMajorParents;
        
        public static ArrayList<Integer> tipClades;
        public static ArrayList<Virus> tipCladeParents;
        
        public static ArrayList<Integer> treeTypes;
        
        public static ArrayList<Double> trunkAntShifts = new ArrayList<Double>();
	
	static final Comparator<Virus> descendantOrder = new Comparator<Virus>() {
		public int compare(Virus v1, Virus v2) {
			Integer descendantsV1 = new Integer(getNumberOfDescendants(v1));
			Integer descendantsV2 = new Integer(getNumberOfDescendants(v2));
			return descendantsV1.compareTo(descendantsV2);
		}
	};
        
        static final Comparator<Virus> birthOrder = new Comparator<Virus>() {
		public int compare(Virus v1, Virus v2) {
			Double birthTimeV1 = v1.getBirth();
			Double birthTimeV2 = v2.getBirth();
			return birthTimeV1.compareTo(birthTimeV2);
		}
	};
		
	// static methods
	public static void add(Virus v) {		
		tips.add(v);
	}
	public static void clear() {
		tips.clear();
	}
	public static List<Virus> getTips() {
		return tips;
	}
	public static Virus getRoot() {
		return root;
	}
	
	// return a random tip that lies between year from and year to
	public static Virus getRandomTipFromTo(double from, double to) {
	
		// fill temporary list
		List<Virus> select = new ArrayList<Virus>();
		for (Virus v : tips) {
			double x = v.getBirth();
			if (x >= from && x < to) {
				select.add(v);
			}
		}
		
		// pull random virus from this list
		Virus rV = null;
		if (select.size() > 0) {	
			int index = Random.nextInt(0,select.size()-1);
			rV = select.get(index);
		}
		return rV;
		
	}
	
	public static int getDemeCount(int d) {
		int count = 0;
		for (Virus v : tips) {
			if (v.getDeme() == d) {
				count++;
			}
		}
		return count;
	}	
		
	// work backwards for each sample filling the children lists
	public static void fillBackward() {
	
		for (Virus child : tips) {
			Virus parent = child.getParent();
			while (parent != null) {
				parent.addChild(child);
				parent.incrementCoverage();
				child = parent;
				parent = child.getParent();
			}
		}
	
	}
	
	public static void dropTips() {
	
		List<Virus> reducedTips = new ArrayList<Virus>();
		for (Virus v : tips) {
			if (Random.nextBoolean(Parameters.treeProportion)) {
				reducedTips.add(v);
			}
		}
		tips = reducedTips;
	
	}

	// marking to by time, not proportional to prevalence
	public static void markTips() {
	
//		for (Virus v : tips) {
//			if (Random.nextBoolean(Parameters.treeProportion)) {
//				while (v.getParent() != null) {
//					v.mark();
//					v = v.getParent();
//				}
//			}
//		}
		
		for (double i = 0; i < Parameters.getDate(); i+=0.1) {
			Virus v = getRandomTipFromTo(i,i+0.1);
			if (v != null) {
				while (v.getParent() != null) {
					v.mark();
					v = v.getParent();
				}
			}
		}
		
	}
		
	// prune tips
	public static void pruneTips() {
	
		List<Virus> reducedTips = new ArrayList<Virus>();
		for (int d = 0; d < Parameters.demeCount; d++) {
			double demeCountTemp = (double) getDemeCount(d);
			double keepProportion = (double) Parameters.tipSamplesPerDeme / (double) getDemeCount(d);
			for (Virus v : tips) {
				if (Random.nextBoolean(keepProportion) && v.getDeme() == d) {
                                        if ((v.getBirth()*365.25) >= Parameters.tipSamplingStartDay & (v.getBirth()*365.25) <= Parameters.tipSamplingEndDay) {
                                            reducedTips.add(v);
                                        }
				}
			}
		}
                samples = new ArrayList<Virus>();
                samples.addAll(tips);
		tips = reducedTips;
                System.out.println("Printing tree with " + tips.size() + " tips");
	
	}
        
        public static void getMRCASeries(ArrayList<Double> dates) {
            
            ArrayList<Virus> activeLineages = new ArrayList<Virus>();
            ArrayList<Virus> sampleQueue = new ArrayList<Virus>();
            sampleQueue.addAll(tips);
            Collections.sort(sampleQueue, birthOrder);
            Collections.reverse(sampleQueue);
            double nextSampleTime = sampleQueue.get(0).getBirth();
            double nextEventTime;
            
            int times = dates.size();
            double currTime; 
            ArrayList<Double> mrcaList = new ArrayList<Double>();
            
            for (int t = times-1; t >= 0; t--) {
                
                currTime = dates.get(t);
                while (nextSampleTime >= currTime) {
                    activeLineages.add(sampleQueue.get(0));
                    sampleQueue.remove(0);
                    if (!sampleQueue.isEmpty())
                        nextSampleTime = sampleQueue.get(0).getBirth();
                    else {
                        nextSampleTime = Double.NEGATIVE_INFINITY;
                    }
                }
                
                Collections.sort(activeLineages, birthOrder); Collections.reverse(activeLineages);
                if (activeLineages.size() > 0) {
                    nextEventTime = activeLineages.get(0).getBirth();
                    while (nextEventTime >= currTime) {
                        Virus pV = activeLineages.get(0).getParent();
                        if (!activeLineages.contains(pV)) {
                            activeLineages.add(pV);
                        }
                        activeLineages.remove(0);
                        Collections.sort(activeLineages, birthOrder); Collections.reverse(activeLineages); //Always need to reorder
                        if (!activeLineages.isEmpty())
                            nextEventTime = activeLineages.get(0).getBirth();
                        else {
                            break;
                        }
                    }
                }
                
                //Get mrca at this time
                int numActive = activeLineages.size();
                double tmrca = 0;
                if (numActive > 1) {
                    for (int i = 0; i < numActive; i++) {
                        for (int j = (i+1); j < numActive; j++) {
                            Virus vA = activeLineages.get(i);
                            Virus vB = activeLineages.get(j);
                            //double dist = vA.distance(vB);
                            double dist = vA.distanceFromTime(vB, currTime);
                            if (dist > tmrca) {
                                tmrca = dist;
                            }
                        }
                    }
                }
                mrcaList.add(tmrca); //should divide by two?
                
            }
            
            Collections.reverse(mrcaList);
            
            String filePath = "logfiles/manual_output_" + Mutantigen.globalSimNum + ".txt";
        	String textToAppend = "About to save the .mrcaSeries file.\n";

        try (FileWriter fw = new FileWriter(filePath, true)) {
            fw.write(textToAppend);
        } catch (IOException e) {
            System.err.println("An error occurred: " + e.getMessage());
        }
        
        
            try {
                String newFileName = "outfiles/out_" + Mutantigen.globalSimNum + ".mrcaSeries";
                File distFile = new File(newFileName);		
                distFile.delete();
                distFile.createNewFile();
                PrintStream stream = new PrintStream(distFile);
                
                for (int t = 0; t < times; t++) {
                    stream.printf("%.4f", dates.get(t));
                    stream.printf("%s", ", ");
                    stream.printf("%.4f", mrcaList.get(t));
                    stream.println();
                }
                  
            } catch(IOException ex) {
                System.out.println("Could not write to file"); 
		System.exit(0);
            }
            
            
        }
        
        public static ArrayList<Integer> getTipTypes() {
            tipTypes = new ArrayList<Integer>();
            for(int t = 0; t < tips.size(); t++) {
                Virus v = tips.get(t);
                Phenotype p = v.getPhenotype();
                int type = p.antigenicType();
                if (!tipTypes.contains(type)) {
                    tipTypes.add(type);
                }
            }
            return tipTypes;
        }
        
        public static ArrayList<Integer> getTreeTypes() {
            treeTypes = new ArrayList<Integer>();
            List<Virus> viruses = postOrderNodes();
            for (Virus v : viruses) {
                Phenotype p = v.getPhenotype();
                int type = p.antigenicType();
                if (!treeTypes.contains(type)) {
                    treeTypes.add(type);
                }
            }
            return treeTypes;
        }
        
        public static void printTrunkAntigenicShifts() {
            
        System.out.println();
        
        String filePath = "logfiles/manual_output_" + Mutantigen.globalSimNum + ".txt";
        	String textToAppend = "About to save the .trunkAntigenicShifts file.\n";

        try (FileWriter fw = new FileWriter(filePath, true)) {
            fw.write(textToAppend);
        } catch (IOException e) {
            System.err.println("An error occurred: " + e.getMessage());
        }
            try {
                String newFileName2 = "outfiles/out_" + Mutantigen.globalSimNum + ".trunkAntigenicShifts";
                File distFile = new File(newFileName2);		
                distFile.delete();
                distFile.createNewFile();
                PrintStream stream = new PrintStream(distFile);
                
                for (Double size : trunkAntShifts) {
                    stream.printf("%.4f", size);
                    stream.println();
                }
                
            } catch(IOException ex) {
                System.out.println("Could not write to file"); 
		System.exit(0);
            }
            
        }
        
        public static ArrayList<Integer> getTipClades() {
            tipClades = new ArrayList<Integer>();   //holds clade of each tip
            tipCladeParents = new ArrayList<Virus>();
            ArrayList<Virus> parents = new ArrayList<Virus>();
            for(int t = 0; t < tips.size(); t++) { //should maybe be for all samples, not tips
                Virus v = tips.get(t);
                Virus trunkParent = getTrunkParent(v);
                if (!tipCladeParents.contains(trunkParent)) {
                    tipCladeParents.add(trunkParent);
                }
                parents.add(trunkParent);
                //tipClades.add(tipCladeParents.indexOf(trunkParent));
            }
            
            Collections.sort((List) tipCladeParents, descendantOrder);
            Collections.reverse(tipCladeParents);
            
            // See how many tips each parent has
            int clades = tipCladeParents.size();
            ArrayList<Integer> counts = new ArrayList<Integer>();
            for(int p = 0; p < clades; p++) {
                counts.add(0);
            }
            for(int t = 0; t < tips.size(); t++) {
                int index = tipCladeParents.indexOf(parents.get(t));
                counts.set(index, counts.get(index)+1);
            }
            
            //ArrayList<Double> birthTimes = new ArrayList<Double>();
            //for (int p = 0; p < clades; p++) {
                //birthTimes.add(tipCladeParents.get(p).getBirth());
            //}
            
            // consolidate clades with less than min tip descendents
            int minCladeSize = 3;
            for (int p = (clades-1); p > 0; p--) {
                int cladeSize = counts.get(p);
                if (cladeSize < minCladeSize) {
                    Virus nextTrunkParent = tipCladeParents.get(p).getParent();
                    while (!tipCladeParents.contains(nextTrunkParent) & nextTrunkParent != null) {
                        nextTrunkParent = nextTrunkParent.getParent();
                    }
                    if (nextTrunkParent == null | !tipCladeParents.contains(nextTrunkParent)) {
                        System.out.println("WARNING: Error traversing trunk!");
                    }
                    int parentIndex = tipCladeParents.indexOf(nextTrunkParent); //p-1, only works if trunk is a single linege
                    counts.remove(p);
                    counts.set(parentIndex, counts.get(parentIndex)+cladeSize);
                    int rndex = parents.indexOf(tipCladeParents.get(p));
                    while (rndex != -1) {
                        parents.set(rndex, tipCladeParents.get(parentIndex));
                        rndex = parents.indexOf(tipCladeParents.get(p));
                    }
                    tipCladeParents.remove(p); 
                }
            }
            
            // need to reassign tip clades after sorting!!!
            for(int t = 0; t < tips.size(); t++) {
                tipClades.add(tipCladeParents.indexOf(parents.get(t)));
            }
                        String filePath = "logfiles/manual_output_" + Mutantigen.globalSimNum + ".txt";
        	String textToAppend = "About to save the .tipClades file.\n";

        try (FileWriter fw = new FileWriter(filePath, true)) {
            fw.write(textToAppend);
        } catch (IOException e) {
            System.err.println("An error occurred: " + e.getMessage());
        } 
            try {
               String newFileName3 = "outfiles/out_" + Mutantigen.globalSimNum + ".tipClades";
                File distFile = new File(newFileName3);		
                distFile.delete();
                distFile.createNewFile();
                PrintStream stream = new PrintStream(distFile);
                
                for (int i = 0; i < tipClades.size(); i++) {
                    stream.printf("%d", tipClades.get(i)+1);
                    if (i < tipClades.size() - 1) {
                        stream.printf("%s", ", ");
                    }
                }
                stream.println();
                
            } catch(IOException ex) {
                System.out.println("Could not write to file"); 
		System.exit(0);
            }
            
            
            return tipClades;
        }
        
        public static Virus getTrunkParent(Virus v) {
            Virus trunkParent;
            if (v.isTrunk()) {
                   trunkParent = v; 
            } else {
                v = v.getParent();
                while (!v.isTrunk()) {
                    v = v.getParent();
                }
                trunkParent = v;
            }
            return trunkParent;
        }
        
        public static void reconstructAntDynamics(ArrayList<Double> dates) {
            
            double timeIncrement = dates.get(dates.size()-1) - dates.get(dates.size()-2);
            
            // set up array for counts
            ArrayList<ArrayList<Integer>> counts = new ArrayList<ArrayList<Integer>>();
            for (int d = 0; d < dates.size(); d++) {
                counts.add(new ArrayList<Integer>());
                for (int j = 0; j <= treeTypes.size(); j++) {
                    counts.get(d).add(0);
                }
            }
            int notInTreeIndex = treeTypes.size();
            
            int slot; int type; int index;
            for (int t = 0; t < samples.size(); t++) { 
                slot = (int) Math.floor(samples.get(t).getBirth() / timeIncrement);    // doesn't work if have a burnin
                if (slot >= dates.size()) {
                    System.out.println("Samples outside of known time indexes");
                    String filePath3 = "logfiles/manual_output_" + Mutantigen.globalSimNum + ".txt";
        		   String textToAppend3 = "Samples outside of known time indexes. \n";

       		try (FileWriter fw = new FileWriter(filePath3, true)) {
           		 fw.write(textToAppend3);
       		} catch (IOException e) {
            		System.err.println("An error occurred: " + e.getMessage());
        		}
                }
                Virus v = samples.get(t);
                Phenotype p = v.getPhenotype();
                type = p.antigenicType(); 
                index = treeTypes.indexOf(type);
                System.out.printf("Index value is: " + index+ "\n");
                if(slot == counts.size()){
                		System.out.printf("code is about to break");
                }
                if (index >= 0) {
                    System.out.printf("last index is:" +  slot + "\n");
                    System.out.println("Length of counts array: " + counts.size());
                    counts.get(slot).set(index, (counts.get(slot).get(index)+1));    
                } else {
                 System.out.printf("last index is:" +  slot + "\n");
                 System.out.println("Length of counts array: " + counts.size());
                    counts.get(slot).set(notInTreeIndex, (counts.get(slot).get(notInTreeIndex)+1));
                }
            }
            
         String filePath = "logfiles/manual_output_" + Mutantigen.globalSimNum + ".txt";
        	String textToAppend = "About to save the .antigenicSamples file.\n";

        try (FileWriter fw = new FileWriter(filePath, true)) {
            fw.write(textToAppend);
        } catch (IOException e) {
            System.err.println("An error occurred: " + e.getMessage());
        }
        
            try {
                String newFileName4 = "outfiles/out_" + Mutantigen.globalSimNum + ".antigenicSamples";
                File seriesFile = new File(newFileName4);		
                seriesFile.delete();
                seriesFile.createNewFile();
                PrintStream stream = new PrintStream(seriesFile);    
                for (int d = 0; d < dates.size(); d++) {
                    stream.printf("%.4f", dates.get(d));
                    stream.printf("%s", ", ");
                    for (int j = 0; j <= treeTypes.size(); j++) {
                        //int currType = majorTypes.get(j);
                        //int currLoc = tipTypes.indexOf(currType);
                        stream.printf("%d", counts.get(d).get(j));
                        if (j < treeTypes.size()) {
                            stream.printf("%s", ", ");
                        }
                    }
                    stream.println();
                }
            } catch(IOException ex) {
                        String filePath1 = "logfiles/manual_output_" + Mutantigen.globalSimNum + ".txt";
        	       String textToAppend1 = "Could not save the .antigenicSamples file.\n";

        try (FileWriter fw = new FileWriter(filePath1, true)) {
            fw.write(textToAppend1);
        } catch (IOException e) {
            System.err.println("An error occurred: " + e.getMessage());
        }
                System.out.println("Could not write to file"); 
		System.exit(0);
            }
            
            
        }
        
        public static void reconstructMajorAntDynamics(ArrayList<Double> dates) {
            
            double timeIncrement = dates.get(dates.size()-1) - dates.get(dates.size()-2);
            
            // set up array for counts
            ArrayList<ArrayList<Integer>> counts = new ArrayList<ArrayList<Integer>>();
            for (int d = 0; d < dates.size(); d++) {
                counts.add(new ArrayList<Integer>());
                for (int j = 0; j < tipTypes.size(); j++) {
                    counts.get(d).add(0);
                }
            }
            
            int slot; int type; int index;
            for (int t = 0; t < tips.size(); t++) { 
                slot = (int) Math.floor(tips.get(t).getBirth() / timeIncrement);    // doesn't work if have a burnin
                if (slot >= dates.size()) {
                    System.out.println("Samples outside of known time indexes");
                }
                Virus v = tips.get(t);
                Phenotype p = v.getPhenotype();
                type = p.antigenicType(); 
                index = tipTypes.indexOf(type);
                if (index == -1) {
                    System.out.println("Type index out of bounds");
                }
                counts.get(slot).set(index, (counts.get(slot).get(index)+1));
            }
            
            majorTypes = new ArrayList<Integer>();
            for (int j = 0; j < tipTypes.size(); j++) {
               double totalCount = 0;
               for (int d = 0; d < dates.size(); d++) {
                   totalCount += counts.get(d).get(j);
               }
               if (totalCount > 20) {
                   majorTypes.add(tipTypes.get(j));
               } 
            }
                        String filePath1 = "logfiles/manual_output_" + Mutantigen.globalSimNum + ".txt";
        	String textToAppend1 = "About to save the .antigenic Samples file for Major Ant dynamics.\n";

        try (FileWriter fw = new FileWriter(filePath1, true)) {
            fw.write(textToAppend1);
        } catch (IOException e) {
            System.err.println("An error occurred: " + e.getMessage());
        } 
            try {
                String newFileName5 = "outfiles/out_" + Mutantigen.globalSimNum + ".antigenicSamples";
                File seriesFile = new File(newFileName5);		
                seriesFile.delete();
                seriesFile.createNewFile();
                PrintStream stream = new PrintStream(seriesFile);    
                for (int d = 0; d < dates.size(); d++) {
                    stream.printf("%.4f", dates.get(d));
                    stream.printf("%s", ", ");
                    for (int j = 0; j < majorTypes.size(); j++) {
                        int currType = majorTypes.get(j);
                        int currLoc = tipTypes.indexOf(currType);
                        stream.printf("%d", counts.get(d).get(currLoc));
                        if (j < majorTypes.size() - 1) {
                            stream.printf("%s", ", ");
                        }
                    }
                    stream.println();
                }
            } catch(IOException ex) {
                        String filePath = "logfiles/manual_output_" + Mutantigen.globalSimNum + ".txt";
        	String textToAppend = "could not save the .antigenicSamples file.\n";

        try (FileWriter fw = new FileWriter(filePath, true)) {
            fw.write(textToAppend);
        } catch (IOException e) {
            System.err.println("An error occurred: " + e.getMessage());
        }
                System.out.println("Could not write to file"); 
		System.exit(0);
            }
            
            
        }
        
        public static void reconstructCladeDynamics(ArrayList<Double> dates) {
            
            double timeIncrement = dates.get(dates.size()-1) - dates.get(dates.size()-2);
            
            // set up array for counts
            ArrayList<ArrayList<Integer>> counts = new ArrayList<ArrayList<Integer>>();
            for (int d = 0; d < dates.size(); d++) {
                counts.add(new ArrayList<Integer>());
                for (int j = 0; j <= tipCladeParents.size(); j++) {
                    counts.get(d).add(0);
                }
            }
            int noCladeIndex = tipCladeParents.size();
            
            int slot; int type; int index;
            for (int t = 0; t < samples.size(); t++) { 
                slot = (int) Math.floor(samples.get(t).getBirth() / timeIncrement);    // doesn't work if have a burnin
                if (slot >= dates.size()) {
                    System.out.println("Samples outside of known time indexes");
                }
                Virus v = samples.get(t);
                Virus trunkParent = getTrunkParent(v);
                index = tipCladeParents.indexOf(trunkParent);
                if (index == -1) {
                    Virus nextTrunkParent = trunkParent.getParent();
                    while (!tipCladeParents.contains(nextTrunkParent) & nextTrunkParent != null) {
                        nextTrunkParent = nextTrunkParent.getParent();
                    }
                    if (nextTrunkParent != null) {
                        index = tipCladeParents.indexOf(nextTrunkParent);
                        counts.get(slot).set(index, (counts.get(slot).get(index)+1));
                    } else {
                        //System.out.println("Trunk parent out of bounds");
                        counts.get(slot).set(noCladeIndex, (counts.get(slot).get(noCladeIndex)+1));
                    }
                } else {
                    counts.get(slot).set(index, (counts.get(slot).get(index)+1));
                }
            }
            
                        String filePath = "logfiles/manual_output_" + Mutantigen.globalSimNum + ".txt";
        	String textToAppend = "About to save the .cladeSamples file.\n";

        try (FileWriter fw = new FileWriter(filePath, true)) {
            fw.write(textToAppend);
        } catch (IOException e) {
            System.err.println("An error occurred: " + e.getMessage());
        }
            
            try {
                String newFileName6 = "outfiles/out_" + Mutantigen.globalSimNum + ".cladeSamples";
                File seriesFile = new File(newFileName6);		
                seriesFile.delete();
                seriesFile.createNewFile();
                PrintStream stream = new PrintStream(seriesFile);    
                for (int d = 0; d < dates.size(); d++) {
                    stream.printf("%.4f", dates.get(d));
                    stream.printf("%s", ", ");
                    for (int j = 0; j < tipCladeParents.size(); j++) {
                        //int currType = majorTypes.get(j);
                        //int currLoc = tipTypes.indexOf(currType);
                        stream.printf("%d", counts.get(d).get(j));
                        if (j < tipCladeParents.size() - 1) {
                            stream.printf("%s", ", ");
                        }
                    }
                    stream.println();
                }
            } catch(IOException ex) {
                System.out.println("Could not write to file"); 
		System.exit(0);
            }
            
            
        }
	
	// returns virus v and all its descendents via a depth-first traversal
	public static List<Virus> postOrderNodes(Virus v) {
		List<Virus> vNodes = new ArrayList<Virus>();
		vNodes.add(v);
		vNodes = postOrderChildren(vNodes);
		return vNodes;
	}
	
	public static List<Virus> postOrderNodes() {
		return postOrderNodes(root);
	}	
	
	// returns virus v and all its descendents via a depth-first traversal
	public static List<Virus> postOrderChildren(List<Virus> vNodes) {
	
		Virus last = vNodes.get(vNodes.size()-1);
	
		for (Virus child : last.getChildren()) {
			vNodes.add(child);
			postOrderChildren(vNodes);
		}
		
		return vNodes;
	
	}


	// Count total descendents of a Virus, working through its children and its children's children
	public static int getNumberOfDescendants(Virus v) {
	
		int numberOfDescendants = v.getNumberOfChildren();
		
		for (Virus child : v.getChildren()) {
			numberOfDescendants += getNumberOfDescendants(child);
		}
		
		return numberOfDescendants;
		
	}
	
	public static int getNumberOfDescendants() {
		return getNumberOfDescendants(root);
	}
		
	// sorts children lists so that first member is child with more descendents than second member
	public static void sortChildrenByDescendants(Virus v) {
		
		List<Virus> children = v.getChildren();
		Collections.sort(children, descendantOrder);
		
		for (Virus child : children) {
			sortChildrenByDescendants(child);
		}
				
	}	
	
	public static void sortChildrenByDescendants() {
		sortChildrenByDescendants(root);
	}
	
	// sets Virus layout based on a postorder traversal
	public static void setLayoutByDescendants() {
	
		List<Virus> vNodes = postOrderNodes();
		
		// set layout of tips based on traversal
		double y = 0;
		for (Virus v : vNodes) {
//			if (tips.contains(v)) {
			if (v.isTip()) {
				v.setLayout(y);
				y++;
			}
		}
		
		// update layout of internal nodes
		Collections.reverse(vNodes);
		for (Virus v : vNodes) {
			if (v.getNumberOfChildren() > 0) {
				double mean = 0;
				for (Virus child : v.getChildren()) {
					mean += child.getLayout();
				}
				mean /= v.getNumberOfChildren();
				v.setLayout(mean);
			}
		}
		
	}	
	
	// looks at a virus and its grandparent, if traits are identical and there is no branching
	// then make virus child rather than grandchild
	// returns v.parent after all is said and done
	public static Virus collapse(Virus v) {
	
		Virus vp = null;
		Virus vgp = null;
		if (v.getParent() != null) {
			vp = v.getParent();
			if (vp.getParent() != null) {
				vgp = vp.getParent();
			}
		}

		if (vp != null && vgp != null) {
	//		if (vp.getNumberOfChildren() == 1 && v.getPhenotype() == vp.getPhenotype() && v.isTrunk() == vp.isTrunk() && v.getDeme() == vp.getDeme()) {
		
			if (vp.getNumberOfChildren() == 1) {
                            // and if tracking viral lineage states
                            Phenotype vPheno = v.getPhenotype();
                            Phenotype vpPheno = v.getPhenotype();
                            int vMutLoad = vPheno.mutLoad();
                            int vpMutLoad = vpPheno.mutLoad();
                            int vAntType = vPheno.antigenicType();
                            int vpAntType = vpPheno.antigenicType();
                            
                            if (vMutLoad == vpMutLoad && vAntType == vpAntType) {

                                    List<Virus> vgpChildren = vgp.getChildren();
                                    int vpIndex =  vgpChildren.indexOf(vp);

                                    if (vpIndex >= 0) {

                                            // replace virus as child of grandparent
                                            vgpChildren.set(vpIndex, v);

                                            // replace grandparent as parent of virus
                                            v.setParent(vgp);

                                            // erase parent
                                            vp = null;

                                    }

                            }
                        }
		}
		
		return v.getParent();

	}
	
	// walks backward using the list of tips, collapsing where possible
	public static void streamline() {
		
		for (Virus v : tips) {
			Virus vp = v;
			while (vp != null) {
				vp = collapse(vp);
			}
		}
		
	}
	
	// rotate the 2d euclidean space using PCA, returning an x-axis with maximum variance
	public static void rotate() {
	
		if (Parameters.phenotypeSpace == "geometric") {
			
			// load a 2d array with phenotypes
			
			List<Virus> virusList = postOrderNodes();
			int n = virusList.size();
			int m = 2;
			
			double[][] input = new double[n][m];
			
			for (int i = 0; i < n; i++) {
				Virus v = virusList.get(i);
				GeometricPhenotype p = (GeometricPhenotype) v.getPhenotype();
				double x = p.getTraitA();
				double y = p.getTraitB();	
				input[i][0] = x;
				input[i][1] = y;				
			}
			
			// project this array
			
			double[][] projected = SimplePCA.project(input);
			
			// reset phenotypes based on projection
			
			for (int i = 0; i < n; i++) {
				Virus v = virusList.get(i);
				GeometricPhenotype p = (GeometricPhenotype) v.getPhenotype();
				double x = projected[i][0];
				double y = projected[i][1];				
				p.setTraitA(x);
				p.setTraitB(y);					
			}

		}	
		
		if (Parameters.phenotypeSpace == "geometric3d") {
			
			// load a 2d array with phenotypes
			
			List<Virus> virusList = postOrderNodes();
			int n = virusList.size();
			int m = 3;
			
			double[][] input = new double[n][m];
			
			for (int i = 0; i < n; i++) {
				Virus v = virusList.get(i);
				GeometricPhenotype3D p = (GeometricPhenotype3D) v.getPhenotype();
				double x = p.getTraitA();
				double y = p.getTraitB();	
				double z = p.getTraitC();	
				input[i][0] = x;
				input[i][1] = y;		
				input[i][2] = z;
			}
			
			// project this array
			
			double[][] projected = SimplePCA.project3D(input);
			
			// reset phenotypes based on projection
			
			for (int i = 0; i < n; i++) {
				Virus v = virusList.get(i);
				GeometricPhenotype3D p = (GeometricPhenotype3D) v.getPhenotype();
				double x = projected[i][0];
				double y = projected[i][1];	
				double z = projected[i][2];	
				p.setTraitA(x);
				p.setTraitB(y);	
				p.setTraitC(z);	
			}

		}			
	
	}
	
	// flips the 2d euclidean space so that first sample is always to the left of the last sample
	public static void flip() {
	
		if (Parameters.phenotypeSpace == "geometric") {

			List<Virus> virusList = postOrderNodes();
			int n = virusList.size();	
			
			// find first and last virus			
			Virus firstVirus = virusList.get(0);
			Virus lastVirus = virusList.get(0);
			double firstDate = firstVirus.getBirth();
			double lastDate = lastVirus.getBirth();
					
			for (Virus v : virusList) {
				if (v.getBirth() < firstDate) {
					firstDate = v.getBirth();
					firstVirus = v;
				}
				if (v.getBirth() > lastDate) {
					lastDate = v.getBirth();
					lastVirus = v;
				}				
			}
			
			// is the x-value of first virus greater than the x-value of last virus?
			// if so, flip
			
			GeometricPhenotype p = (GeometricPhenotype) firstVirus.getPhenotype();
			double firstX = p.getTraitA();
			p = (GeometricPhenotype) lastVirus.getPhenotype();
			double lastX = p.getTraitA();		
			
			if (firstX > lastX) {
			
				// I think that postOrderNodes() has replicates in it, need to go through some hoops because of this
				double[] input = new double[n];
			
				for (int i = 0; i < n; i++) {
					Virus v = virusList.get(i);
					p = (GeometricPhenotype) v.getPhenotype();		
					input[i] = p.getTraitA();;
				}
				
				for (int i = 0; i < n; i++) {
					Virus v = virusList.get(i);
					p = (GeometricPhenotype) v.getPhenotype();
					double x = -1*input[i];		
					p.setTraitA(x);
				}				
			
			}
			
		}
		
		if (Parameters.phenotypeSpace == "geometric3d") {

			List<Virus> virusList = postOrderNodes();
			int n = virusList.size();	
			
			// find first and last virus			
			Virus firstVirus = virusList.get(0);
			Virus lastVirus = virusList.get(0);
			double firstDate = firstVirus.getBirth();
			double lastDate = lastVirus.getBirth();
					
			for (Virus v : virusList) {
				if (v.getBirth() < firstDate) {
					firstDate = v.getBirth();
					firstVirus = v;
				}
				if (v.getBirth() > lastDate) {
					lastDate = v.getBirth();
					lastVirus = v;
				}				
			}
			
			// is the x-value of first virus greater than the x-value of last virus?
			// if so, flip
			
			GeometricPhenotype3D p = (GeometricPhenotype3D) firstVirus.getPhenotype();
			double firstX = p.getTraitA();
			p = (GeometricPhenotype3D) lastVirus.getPhenotype();
			double lastX = p.getTraitA();		
			
			if (firstX > lastX) {
			
				// I think that postOrderNodes() has replicates in it, need to go through some hoops because of this
				double[] input = new double[n];
			
				for (int i = 0; i < n; i++) {
					Virus v = virusList.get(i);
					p = (GeometricPhenotype3D) v.getPhenotype();		
					input[i] = p.getTraitA();;
				}
				
				for (int i = 0; i < n; i++) {
					Virus v = virusList.get(i);
					p = (GeometricPhenotype3D) v.getPhenotype();
					double x = -1*input[i];		
					p.setTraitA(x);
				}				
			
			}
			
		}		
	
	}
	
	// walks through list of nodes and update min and max ranges appropriately
	public static void updateRange() {
	
		xMin = 0.0;
		xMax = 0.0;
		yMin = 0.0;
		yMax = 0.0;
		zMin = 0.0;
		zMax = 0.0;		
	
		if (Parameters.phenotypeSpace == "geometric") {
			for (Virus v : postOrderNodes()) {
			
				GeometricPhenotype p = (GeometricPhenotype) v.getPhenotype();
				double x = p.getTraitA();
				double y = p.getTraitB();
				if (xMin > x) { xMin = x; }
				if (xMax < x) { xMax = x; }
				if (yMin > y) { yMin = y; }
				if (yMax < y) { yMax = y; }	
			
			}
		}
		
		if (Parameters.phenotypeSpace == "geometric3d") {
			for (Virus v : postOrderNodes()) {
			
				GeometricPhenotype3D p = (GeometricPhenotype3D) v.getPhenotype();
				double x = p.getTraitA();
				double y = p.getTraitB();
				double z = p.getTraitC();				
				if (xMin > x) { xMin = x; }
				if (xMax < x) { xMax = x; }
				if (yMin > y) { yMin = y; }
				if (yMax < y) { yMax = y; }	
				if (zMin > z) { zMin = z; }
				if (zMax < z) { zMax = z; }					
			
			}
		}		
		
		xMin = Math.floor(xMin) - 10;
		xMax = Math.ceil(xMax) + 10;
		yMin = Math.floor(yMin) - 10;
		yMax = Math.ceil(yMax) + 10;
		zMin = Math.floor(zMin) - 10;
		zMax = Math.ceil(zMax) + 10;		
	
	}

	public static void printRange() {
		            String filePath = "logfiles/manual_output_" + Mutantigen.globalSimNum + ".txt";
        	String textToAppend = "About to save the .range file.\n";

        try (FileWriter fw = new FileWriter(filePath, true)) {
            fw.write(textToAppend);
        } catch (IOException e) {
            System.err.println("An error occurred: " + e.getMessage());
        }
		try {
		        String newFileName7 = "outfiles/out_" + Mutantigen.globalSimNum + ".range";
			File rangeFile = new File(newFileName7);
			rangeFile.delete();
			rangeFile.createNewFile();
			PrintStream rangeStream = new PrintStream(rangeFile);
			rangeStream.printf("%.4f,%.4f,%.4f,%.4f,%.4f,%.4f\n", xMin, xMax, yMin, yMax, zMin, zMax);
			rangeStream.close();
		} catch(IOException ex) {
			System.out.println("Could not write to file"); 
			System.exit(0);
		}
		
	}
	
	public static void printTips() {
		
		            String filePath = "logfiles/manual_output_" + Mutantigen.globalSimNum + ".txt";
        	String textToAppend = "About to save the .tips file.\n";

        try (FileWriter fw = new FileWriter(filePath, true)) {
            fw.write(textToAppend);
        } catch (IOException e) {
            System.err.println("An error occurred: " + e.getMessage());
        }
		try {
			String newFileName8 = "outfiles/out_" + Mutantigen.globalSimNum + ".tips";
			File tipFile = new File(newFileName8);
			tipFile.delete();
			tipFile.createNewFile();
			PrintStream tipStream = new PrintStream(tipFile);
			tipStream.printf("\"%s\",\"%s\",\"%s\",\"%s\",\"%s\",\"%s\",\"%s\",\"%s\",\"%s\"\n", "name", "year", "trunk", "tip", "mark", "location", "layout", "ag1", "ag2");
			for (int i = 0; i < tips.size(); i++) {
				Virus v = tips.get(i);			
				tipStream.printf("\"%s\",%.4f,%d,%d,%d,%d,%.4f,%s\n", v, v.getBirth(), v.isTrunk()?1:0, v.isTip()?1:0, v.isMarked()?1:0, v.getDeme(), v.getLayout(), v.getPhenotype());
			}
			tipStream.close();
		} catch(IOException ex) {
			System.out.println("Could not write to file"); 
			System.exit(0);
		}
		
	}
	
	public static void printBranches() {
		            String filePath = "logfiles/manual_output_" + Mutantigen.globalSimNum + ".txt";
        	String textToAppend = "About to save the .branches file.\n";

        try (FileWriter fw = new FileWriter(filePath, true)) {
            fw.write(textToAppend);
        } catch (IOException e) {
            System.err.println("An error occurred: " + e.getMessage());
        }
		try {
			String newFileName9 = "outfiles/out_" + Mutantigen.globalSimNum + ".branches";
			File branchFile = new File(newFileName9);
			branchFile.delete();
			branchFile.createNewFile();
			PrintStream branchStream = new PrintStream(branchFile);
			for (Virus v : postOrderNodes()) {
				if (v.getParent() != null) {
					Virus vp = v.getParent();
					branchStream.printf("{\"%s\",%.4f,%d,%d,%d,%d,%.4f,%s}\t", v, v.getBirth(), v.isTrunk()?1:0, v.isTip()?1:0, v.isMarked()?1:0, v.getDeme(), v.getLayout(), v.getPhenotype());
					branchStream.printf("{\"%s\",%.4f,%d,%d,%d,%d,%.4f,%s}\t", vp, vp.getBirth(), vp.isTrunk()?1:0, vp.isTip()?1:0, v.isMarked()?1:0, vp.getDeme(), vp.getLayout(), vp.getPhenotype());
					branchStream.printf("%d\n", vp.getCoverage());
				}
			}
			branchStream.close();
		} catch(IOException ex) {
			System.out.println("Could not write to file"); 
			System.exit(0);
		}
		
	}
	
	// assess node in building Newick string
	public static Virus assessNode(Virus v, List<Virus> visited, PrintStream treeStream) {
	
		Virus returnVirus = null;
		boolean printHeight = false;
	
		// if virus has multiple children, return first child that has not been visited
		if (v.getNumberOfChildren() > 1) {
			boolean childrenVisited = true;
			for (int i = 0; i < v.getNumberOfChildren(); i++) {
				Virus vc = v.getChildren().get(i);
				if (!visited.contains(vc)) {
					if (i == 0) {
						treeStream.print("(");
					}
					else {
						treeStream.print(",");
					}
					childrenVisited = false;
					returnVirus = vc;
					break;
				}
			}
			// failure, all children visited, return to parent
			if (childrenVisited) {
				treeStream.print(")");	
				printHeight = true;
				returnVirus = v.getParent();
			}
		}
		
		// if tip is encountered, print tip, return to parent
		if (v.getNumberOfChildren() == 0) {
			//String nameOut = v.toString() + Double.toString(v.getBirth());
                        treeStream.print(v.toString());
			printHeight = true;
			returnVirus = v.getParent();		
		}			
		
		// walk down (or up) branches 
		if (v.getNumberOfChildren() == 1) {
			Virus vc = v.getChildren().get(0);
			if (!visited.contains(vc)) {
				returnVirus = vc;
			}
			else {
				returnVirus = v.getParent();
			}
		}
		
		// find height, walk back until a parent with a split occurs
		if (printHeight && v.getParent() != null) {

			//treeStream.printf("[&antigenic={%s}]", v.getPhenotype() );
		
			Virus vp = v.getParent();
			while (vp.getNumberOfChildren() == 1 && vp.getParent() != null) {
				vp = vp.getParent();
			}
			double height = v.getBirth() - vp.getBirth();
			treeStream.printf(":%.4f", height);	

		}
		
		return returnVirus;
	
	}
        
        	// assess node in building Newick string
	public static Virus assessNodeSimMap(Virus v, List<Virus> visited, PrintStream treeStream) {
	
		Virus returnVirus = null;
		boolean printHeight = false;
	
		// if virus has multiple children, return first child that has not been visited
                
                if (v.getNumberOfChildren() > 2) {
                    System.out.println("Hit non-bifurcating node with " + v.getNumberOfChildren() + " chilren!");
                }
                
		if (v.getNumberOfChildren() > 1) {
			boolean childrenVisited = true;
			for (int i = 0; i < v.getNumberOfChildren(); i++) {
				Virus vc = v.getChildren().get(i);
				if (!visited.contains(vc)) {
					if (i == 0) {
						treeStream.print("(");
					}
					else {
						treeStream.print(",");
					}
					childrenVisited = false;
					returnVirus = vc;
					break;
				}
			}
			// failure, all children visited, return to parent
			if (childrenVisited) {
				treeStream.print(")");	
				printHeight = true;
				returnVirus = v.getParent();
			}
		}
		
		// if tip is encountered, print tip, return to parent
		if (v.getNumberOfChildren() == 0) {
			treeStream.print(v.toString());
			printHeight = true;
			returnVirus = v.getParent();		
		}			
		
		// walk down (or up) branches 
		if (v.getNumberOfChildren() == 1) {
			Virus vc = v.getChildren().get(0);
			if (!visited.contains(vc)) {
				returnVirus = vc;
			}
			else {
				returnVirus = v.getParent();
			}
		}
		
		// find height, walk back until a parent with a split occurs
		if (printHeight && v.getParent() != null) {

			//treeStream.printf("[&antigenic={%s}]", v.getPhenotype() );
                        String lineString = "{";
                    
			Virus vp = v.getParent();
                        double segTime = v.getBirth() - vp.getBirth();
                        double heightCheck = segTime;
                        Phenotype p = v.getPhenotype();
                        int state = p.mutLoad();
                        lineString = lineString + Integer.toString(state) + "," + Double.toString(segTime);
                        Virus vc;
			while (vp.getNumberOfChildren() == 1 && vp.getParent() != null) {
                            
                            vc = vp;
                            vp = vp.getParent();
                            
                            lineString = lineString + ":";
                            segTime = vc.getBirth() - vp.getBirth();
                            heightCheck += segTime;
                            p = vc.getPhenotype();
                            state = p.mutLoad();
                            lineString = lineString + Integer.toString(state) + "," + Double.toString(segTime);
                            
			}
                        
                        lineString = lineString + "}";
                        //System.out.println(lineString);
                        
			double height = v.getBirth() - vp.getBirth();
                        if (height <= 0) {
                            System.out.println("Hit zero or negative branch lengths!");
                        }
                        //System.out.println("Expected height = " + Double.toString(height) + "; SimMap height = " + Double.toString(heightCheck));
			//treeStream.printf(":%.4f", height);
                        treeStream.printf(":%s", lineString);

		}
		
		return returnVirus;
	
	}
        
        public static Virus assessNodeAntigenic(Virus v, List<Virus> visited, PrintStream treeStream) {
	
		Virus returnVirus = null;
		boolean printHeight = false;
	
		// if virus has multiple children, return first child that has not been visited
                
                if (v.getNumberOfChildren() > 2) {
                    System.out.println("Hit non-bifurcating node with " + v.getNumberOfChildren() + " chilren!");
                }
                
		if (v.getNumberOfChildren() > 1) {
			boolean childrenVisited = true;
			for (int i = 0; i < v.getNumberOfChildren(); i++) {
				Virus vc = v.getChildren().get(i);
				if (!visited.contains(vc)) {
					if (i == 0) {
						treeStream.print("(");
					}
					else {
						treeStream.print(",");
					}
					childrenVisited = false;
					returnVirus = vc;
					break;
				}
			}
			// failure, all children visited, return to parent
			if (childrenVisited) {
				treeStream.print(")");	
				printHeight = true;
				returnVirus = v.getParent();
			}
		}
		
		// if tip is encountered, print tip, return to parent
		if (v.getNumberOfChildren() == 0) {
			treeStream.print(v.toString());
			printHeight = true;
			returnVirus = v.getParent();		
		}			
		
		// walk down (or up) branches 
		if (v.getNumberOfChildren() == 1) {
			Virus vc = v.getChildren().get(0);
			if (!visited.contains(vc)) {
				returnVirus = vc;
			}
			else {
				returnVirus = v.getParent();
			}
		}
		
		// find height, walk back until a parent with a split occurs
		if (printHeight && v.getParent() != null) {

			//treeStream.printf("[&antigenic={%s}]", v.getPhenotype() );
                        String lineString = "{";
                    
			Virus vp = v.getParent();
                        double segTime = v.getBirth() - vp.getBirth();
                        double heightCheck = segTime;
                        Phenotype p = v.getPhenotype();
                        int type = p.antigenicType();
                        if (!treeTypes.contains(type)) {
                            System.out.println("Antigenic type not found in tree");
                        }
                        int state = treeTypes.indexOf(type) + 1;
                        lineString = lineString + Integer.toString(state) + "," + Double.toString(segTime);
                        Virus vc;
			while (vp.getNumberOfChildren() == 1 && vp.getParent() != null) {
                            
                            vc = vp;
                            vp = vp.getParent();
                            
                            lineString = lineString + ":";
                            segTime = vc.getBirth() - vp.getBirth();
                            heightCheck += segTime;
                            p = vc.getPhenotype();
                            type = p.antigenicType();
                            if (!treeTypes.contains(type)) {
                                System.out.println("Antigenic type not found in tree");
                            }
                            state = treeTypes.indexOf(type) + 1;
                            lineString = lineString + Integer.toString(state) + "," + Double.toString(segTime);
                            
			}
                        
                        lineString = lineString + "}";
                        //System.out.println(lineString);
                        
			double height = v.getBirth() - vp.getBirth();
                        if (height <= 0) {
                            System.out.println("Hit zero or negative branch lengths!");
                        }
                        //System.out.println("Expected height = " + Double.toString(height) + "; SimMap height = " + Double.toString(heightCheck));
			//treeStream.printf(":%.4f", height);
                        treeStream.printf(":%s", lineString);

		}
		
		return returnVirus;
	
	}
        
        public static Virus assessNodeAntigenicTransitions(Virus v, List<Virus> visited, PrintStream treeStream) {
	
		Virus returnVirus = null;
		boolean printHeight = false;
	
		// if virus has multiple children, return first child that has not been visited
                
                if (v.getNumberOfChildren() > 2) {
                    System.out.println("Hit non-bifurcating node with " + v.getNumberOfChildren() + " chilren!");
                }
                
		if (v.getNumberOfChildren() > 1) {
			boolean childrenVisited = true;
			for (int i = 0; i < v.getNumberOfChildren(); i++) {
				Virus vc = v.getChildren().get(i);
				if (!visited.contains(vc)) {
					if (i == 0) {
						treeStream.print("(");
					}
					else {
						treeStream.print(",");
					}
					childrenVisited = false;
					returnVirus = vc;
					break;
				}
			}
			// failure, all children visited, return to parent
			if (childrenVisited) {
				treeStream.print(")");	
				printHeight = true;
				returnVirus = v.getParent();
			}
		}
		
		// if tip is encountered, print tip, return to parent
		if (v.getNumberOfChildren() == 0) {
			treeStream.print(v.toString());
			printHeight = true;
			returnVirus = v.getParent();		
		}			
		
		// walk down (or up) branches 
		if (v.getNumberOfChildren() == 1) {
			Virus vc = v.getChildren().get(0);
			if (!visited.contains(vc)) {
				returnVirus = vc;
			}
			else {
				returnVirus = v.getParent();
			}
		}
		
		// find height, walk back until a parent with a split occurs
		if (printHeight && v.getParent() != null) {

			//treeStream.printf("[&antigenic={%s}]", v.getPhenotype() );
                        String lineString = "{";
                    
			Virus vp = v.getParent();
                        double segTime = v.getBirth() - vp.getBirth();
                        double heightCheck = segTime;
                        Phenotype p = v.getPhenotype();
                        int type = p.antigenicType();
                        if (!treeTypes.contains(type)) {
                            System.out.println("Antigenic type not found in tree");
                        }
                        int state = treeTypes.indexOf(type) + 1;
                        lineString = lineString + Integer.toString(state) + "," + Double.toString(segTime);
                        Virus vc;
                        int currType = type;
			while (vp.getNumberOfChildren() == 1 && vp.getParent() != null) {
                            
                            vc = vp;
                            vp = vp.getParent();
                            
                            lineString = lineString + ":";
                            segTime = vc.getBirth() - vp.getBirth();
                            heightCheck += segTime;
                            p = vc.getPhenotype();
                            type = p.antigenicType();
                            if (type != currType & vc.isTrunk()) {
                                double size = AntigenicTree.computeDistance(currType, type);
                                trunkAntShifts.add(size);
                                currType = type;    
                            }
                            if (!treeTypes.contains(type)) {
                                System.out.println("Antigenic type not found in tree");
                            }
                            state = treeTypes.indexOf(type) + 1;
                            lineString = lineString + Integer.toString(state) + "," + Double.toString(segTime);
                            
			}
                        
                        if (vp.getParent() != null) {
                            p = vp.getPhenotype();
                            type = p.antigenicType();
                            if (type != currType & vp.isTrunk()) {
                                double size = AntigenicTree.computeDistance(currType, type);
                                trunkAntShifts.add(size);
                                //System.out.println("Antigenic transition at internal node!");
                                //currType = type;    
                            }
                        }
                        
                        lineString = lineString + "}";
                        //System.out.println(lineString);
                        
			double height = v.getBirth() - vp.getBirth();
                        if (height <= 0) {
                            System.out.println("Hit zero or negative branch lengths!");
                        }
                        //System.out.println("Expected height = " + Double.toString(height) + "; SimMap height = " + Double.toString(heightCheck));
			//treeStream.printf(":%.4f", height);
                        treeStream.printf(":%s", lineString);

		}
		
		return returnVirus;
	
	}
        
        public static Virus assessNodeMajorAntigenic(Virus v, List<Virus> visited, PrintStream treeStream) {
	
		Virus returnVirus = null;
		boolean printHeight = false;
	
		// if virus has multiple children, return first child that has not been visited
                
                if (v.getNumberOfChildren() > 2) {
                    System.out.println("Hit non-bifurcating node with " + v.getNumberOfChildren() + " chilren!");
                }
                
		if (v.getNumberOfChildren() > 1) {
			boolean childrenVisited = true;
			for (int i = 0; i < v.getNumberOfChildren(); i++) {
				Virus vc = v.getChildren().get(i);
				if (!visited.contains(vc)) {
					if (i == 0) {
						treeStream.print("(");
					}
					else {
						treeStream.print(",");
					}
					childrenVisited = false;
					returnVirus = vc;
					break;
				}
			}
			// failure, all children visited, return to parent
			if (childrenVisited) {
				treeStream.print(")");	
				printHeight = true;
				returnVirus = v.getParent();
			}
		}
		
		// if tip is encountered, print tip, return to parent
		if (v.getNumberOfChildren() == 0) {
			treeStream.print(v.toString());
			printHeight = true;
			returnVirus = v.getParent();		
		}			
		
		// walk down (or up) branches 
		if (v.getNumberOfChildren() == 1) {
			Virus vc = v.getChildren().get(0);
			if (!visited.contains(vc)) {
				returnVirus = vc;
			}
			else {
				returnVirus = v.getParent();
			}
		}
		
		// find height, walk back until a parent with a split occurs
		if (printHeight && v.getParent() != null) {

			//treeStream.printf("[&antigenic={%s}]", v.getPhenotype() );
                        String lineString = "{";
                    
			Virus vp = v.getParent();
                        double segTime = v.getBirth() - vp.getBirth();
                        double heightCheck = segTime;
                        Phenotype p = v.getPhenotype();
                        int type = p.antigenicType();
                        if (!majorTypes.contains(type)) {
                            type = firstMajorParents.get(type);
                        }
                        int state = majorTypes.indexOf(type) + 1;
                        lineString = lineString + Integer.toString(state) + "," + Double.toString(segTime);
                        Virus vc;
			while (vp.getNumberOfChildren() == 1 && vp.getParent() != null) {
                            
                            vc = vp;
                            vp = vp.getParent();
                            
                            lineString = lineString + ":";
                            segTime = vc.getBirth() - vp.getBirth();
                            heightCheck += segTime;
                            p = vc.getPhenotype();
                            type = p.antigenicType();
                            if (!majorTypes.contains(type)) {
                                type = firstMajorParents.get(type);
                            }
                            state = majorTypes.indexOf(type) + 1;
                            lineString = lineString + Integer.toString(state) + "," + Double.toString(segTime);
                            
			}
                        
                        lineString = lineString + "}";
                        //System.out.println(lineString);
                        
			double height = v.getBirth() - vp.getBirth();
                        if (height <= 0) {
                            System.out.println("Hit zero or negative branch lengths!");
                        }
                        //System.out.println("Expected height = " + Double.toString(height) + "; SimMap height = " + Double.toString(heightCheck));
			//treeStream.printf(":%.4f", height);
                        treeStream.printf(":%s", lineString);

		}
		
		return returnVirus;
	
	}
        
        public static Virus assessNodeClade(Virus v, List<Virus> visited, PrintStream treeStream) {
	
		Virus returnVirus = null;
		boolean printHeight = false;
	
		// if virus has multiple children, return first child that has not been visited
                
                if (v.getNumberOfChildren() > 2) {
                    System.out.println("Hit non-bifurcating node with " + v.getNumberOfChildren() + " chilren!");
                }
                
		if (v.getNumberOfChildren() > 1) {
			boolean childrenVisited = true;
			for (int i = 0; i < v.getNumberOfChildren(); i++) {
				Virus vc = v.getChildren().get(i);
				if (!visited.contains(vc)) {
					if (i == 0) {
						treeStream.print("(");
					}
					else {
						treeStream.print(",");
					}
					childrenVisited = false;
					returnVirus = vc;
					break;
				}
			}
			// failure, all children visited, return to parent
			if (childrenVisited) {
				treeStream.print(")");	
				printHeight = true;
				returnVirus = v.getParent();
			}
		}
		
		// if tip is encountered, print tip, return to parent
		if (v.getNumberOfChildren() == 0) {
			treeStream.print(v.toString());
			printHeight = true;
			returnVirus = v.getParent();		
		}			
		
		// walk down (or up) branches 
		if (v.getNumberOfChildren() == 1) {
			Virus vc = v.getChildren().get(0);
			if (!visited.contains(vc)) {
				returnVirus = vc;
			}
			else {
				returnVirus = v.getParent();
			}
		}
		
		// find height, walk back until a parent with a split occurs
		if (printHeight && v.getParent() != null) {

			//treeStream.printf("[&antigenic={%s}]", v.getPhenotype() );
                        String lineString = "{";
                    
			Virus vp = v.getParent();
                        double segTime = v.getBirth() - vp.getBirth();
                        double heightCheck = segTime;
                        
                        Virus trunkParent = getTrunkParent(v);
                        int index = tipCladeParents.indexOf(trunkParent);
                        if (index == -1) {
                            Virus nextTrunkParent = trunkParent.getParent();
                            while (!tipCladeParents.contains(nextTrunkParent) & nextTrunkParent != null) {
                                nextTrunkParent = nextTrunkParent.getParent();
                            }
                            if (nextTrunkParent != null) {
                                index = tipCladeParents.indexOf(nextTrunkParent);
                            } else {
                                //System.out.println("Trunk parent out of bounds");
                            }
                        }
                        int state = index + 1;
                        lineString = lineString + Integer.toString(state) + "," + Double.toString(segTime);
                        Virus vc;
			while (vp.getNumberOfChildren() == 1 && vp.getParent() != null) {
                            
                            vc = vp;
                            vp = vp.getParent();
                            
                            lineString = lineString + ":";
                            segTime = vc.getBirth() - vp.getBirth();
                            heightCheck += segTime;
                            
                            trunkParent = getTrunkParent(vc);
                            index = tipCladeParents.indexOf(trunkParent);
                            if (index == -1) {
                                Virus nextTrunkParent = trunkParent.getParent();
                                while (!tipCladeParents.contains(nextTrunkParent) & nextTrunkParent != null) {
                                    nextTrunkParent = nextTrunkParent.getParent();
                                }
                                if (nextTrunkParent != null) {
                                    index = tipCladeParents.indexOf(nextTrunkParent);
                                } else {
                                    //System.out.println("Trunk parent out of bounds");
                                }
                            }
                            state = index + 1;
                            lineString = lineString + Integer.toString(state) + "," + Double.toString(segTime);
                            
			}
                        
                        lineString = lineString + "}";
                        //System.out.println(lineString);
                        
			double height = v.getBirth() - vp.getBirth();
                        if (height <= 0) {
                            System.out.println("Hit zero or negative branch lengths!");
                        }
                        //System.out.println("Expected height = " + Double.toString(height) + "; SimMap height = " + Double.toString(heightCheck));
			//treeStream.printf(":%.4f", height);
                        treeStream.printf(":%s", lineString);

		}
		
		return returnVirus;
	
	}
	
	public static void printNewick() {
	
	            String filePath = "logfiles/manual_output_" + Mutantigen.globalSimNum + ".txt";
        	String textToAppend = "About to save the .trees file.\n";

        try (FileWriter fw = new FileWriter(filePath, true)) {
            fw.write(textToAppend);
        } catch (IOException e) {
            System.err.println("An error occurred: " + e.getMessage());
        }
        
		try {
			String newFileName10 = "outfiles/out_" + Mutantigen.globalSimNum + ".trees";
			File treeFile = new File(newFileName10);
			treeFile.delete();
			treeFile.createNewFile();
			PrintStream treeStream = new PrintStream(treeFile);
			
			List<Virus> visited = new ArrayList<Virus>();
				
			// start at root
			Virus v = root;
			visited.add(v);			
			
			while (v != null) {
				
				v = assessNode(v, visited, treeStream);		
				visited.add(v);
							
			}
			treeStream.println();
			
			treeStream.close();
		} catch(IOException ex) {
			System.out.println("Could not write to file"); 
			System.exit(0);
		}
	
	}
        
        public static void printSimMapLoad() {
	
	            String filePath = "logfiles/manual_output_" + Mutantigen.globalSimNum + ".txt";
        	String textToAppend = "About to save the .simmapLoad file.\n";

        try (FileWriter fw = new FileWriter(filePath, true)) {
            fw.write(textToAppend);
        } catch (IOException e) {
            System.err.println("An error occurred: " + e.getMessage());
        }
	
		try {
			String newFileName11 = "outfiles/out_" + Mutantigen.globalSimNum + ".simmapLoad";
			File treeFile = new File(newFileName11);
			treeFile.delete();
			treeFile.createNewFile();
			PrintStream treeStream = new PrintStream(treeFile);
			
			List<Virus> visited = new ArrayList<Virus>();
				
			// start at root
			Virus v = root;
			visited.add(v);			
			
			while (v != null) {
				
				v = assessNodeSimMap(v, visited, treeStream);		
				visited.add(v);
							
			}
			treeStream.println();
			
			treeStream.close();
		} catch(IOException ex) {
			System.out.println("Could not write to file"); 
			System.exit(0);
		}
	
	}
        
        public static void printSimMapAntigenic() {
	
		            String filePath = "logfiles/manual_output_" + Mutantigen.globalSimNum + ".txt";
        	String textToAppend = "About to save the .simmapAntigenic file.\n";

        try (FileWriter fw = new FileWriter(filePath, true)) {
            fw.write(textToAppend);
        } catch (IOException e) {
            System.err.println("An error occurred: " + e.getMessage());
        }
	
		try {
			String newFileName12 = "outfiles/out_" + Mutantigen.globalSimNum + ".out.simmapAntigenic";
			File treeFile = new File(newFileName12);
			treeFile.delete();
			treeFile.createNewFile();
			PrintStream treeStream = new PrintStream(treeFile);
                        
			List<Virus> visited = new ArrayList<Virus>();
				
			// start at root
			Virus v = root;
			visited.add(v);			
			
			while (v != null) {
				
				v = assessNodeAntigenicTransitions(v, visited, treeStream);		
				visited.add(v);
							
			}
			treeStream.println();
			
			treeStream.close();
		} catch(IOException ex) {
			System.out.println("Could not write to file"); 
			System.exit(0);
		}
	
	}
        
        public static void printSimMapMajorAntigenic() {
	
	            String filePath = "logfiles/manual_output_" + Mutantigen.globalSimNum + ".txt";
        	String textToAppend = "About to save the .simmapAntigenic file within the major antigenic function .\n";

        try (FileWriter fw = new FileWriter(filePath, true)) {
            fw.write(textToAppend);
        } catch (IOException e) {
            System.err.println("An error occurred: " + e.getMessage());
        }
        
		try {
			String newFileName13 = "outfiles/out_" + Mutantigen.globalSimNum + ".simmapAntigenic";
			File treeFile = new File(newFileName13);
			treeFile.delete();
			treeFile.createNewFile();
			PrintStream treeStream = new PrintStream(treeFile);
			
                        setFirstMajorParents();
                        
			List<Virus> visited = new ArrayList<Virus>();
				
			// start at root
			Virus v = root;
			visited.add(v);			
			
			while (v != null) {
				
				v = assessNodeAntigenic(v, visited, treeStream);		
				visited.add(v);
							
			}
			treeStream.println();
			
			treeStream.close();
		} catch(IOException ex) {
			System.out.println("Could not write to file"); 
			System.exit(0);
		}
	
	}
        
        public static void printSimMapClades() {
	
	            String filePath = "logfiles/manual_output_" + Mutantigen.globalSimNum + ".txt";
        	String textToAppend = "About to save the .simmapClades file.\n";

        try (FileWriter fw = new FileWriter(filePath, true)) {
            fw.write(textToAppend);
        } catch (IOException e) {
            System.err.println("An error occurred: " + e.getMessage());
        }
        
		try {
			String newFileName14 = "outfiles/out_" + Mutantigen.globalSimNum + ".simmapClades";
			File treeFile = new File(newFileName14);
			treeFile.delete();
			treeFile.createNewFile();
			PrintStream treeStream = new PrintStream(treeFile);
			
                        //setFirstMajorParents();
                        
			List<Virus> visited = new ArrayList<Virus>();
				
			// start at root
			Virus v = root;
			visited.add(v);			
			
			while (v != null) {
				
				v = assessNodeClade(v, visited, treeStream);		
				visited.add(v);
							
			}
			treeStream.println();
			
			treeStream.close();
		} catch(IOException ex) {
			System.out.println("Could not write to file"); 
			System.exit(0);
		}
	
	}
        
        public static void setFirstMajorParents() {
            
            firstMajorParents = new ArrayList<Integer>();
            int allTypes = AntigenicTree.nodeCount();
            for (int i = 0; i < allTypes; i++) {
                if (majorTypes.contains(i)) {
                    firstMajorParents.add(i);
                } else {
                    ArrayList<Integer> ancestry = AntigenicTree.getAncestry(i);
                    if (ancestry.size() > 1) {
                        int nextIndex = 1;
                        int nextParent = ancestry.get(nextIndex);
                        while (!majorTypes.contains(nextParent)) {
                            nextIndex++;
                            if (nextIndex >= ancestry.size()) {
                                firstMajorParents.add(0);
                            } else {
                                nextParent = ancestry.get(nextIndex);
                            }
                        }
                        firstMajorParents.add(nextParent);
                    } else {
                        firstMajorParents.add(0);
                    }
                }
            }
        } 
	
	public static int sideBranchMutations() {
		int count = 0;
		for (Virus v : postOrderNodes()) {
			if (v.getParent() != null && v.getBirth() < Parameters.getDate() - Parameters.yearsFromMK) {
				Virus vp = v.getParent();
				if (!v.isTrunk() && !vp.isTrunk() && v.getPhenotype() != vp.getPhenotype()) {
					count++;
				}
			}
		}
		return count;
	}
        
	
	public static double sideBranchOpportunity() {
		double time = 0;
		for (Virus v : postOrderNodes()) {
			if (v.getParent() != null && v.getBirth() < Parameters.getDate() - Parameters.yearsFromMK) {
				Virus vp = v.getParent();
				if (!v.isTrunk() && !vp.isTrunk()) {
					time += v.getBirth() - vp.getBirth();
				}
			}
		}
		return time;
	}	
	
	public static int trunkMutations() {
		int count = 0;
		for (Virus v : postOrderNodes()) {
			if (v.getParent() != null && v.getBirth() < Parameters.getDate() - Parameters.yearsFromMK) {
				Virus vp = v.getParent();
				if (v.isTrunk() && vp.isTrunk() && v.getPhenotype() != vp.getPhenotype()) {
					count++;
				}
			}
		}
		return count;
	}	
	
	public static double trunkOpportunity() {
		double time = 0;
		for (Virus v : postOrderNodes()) {
			if (v.getParent() != null && v.getBirth() < Parameters.getDate() - Parameters.yearsFromMK) {
				Virus vp = v.getParent();
				if (v.isTrunk() && vp.isTrunk()) {
					time += v.getBirth() - vp.getBirth();
				}
			}
		}
		return time;
	}		
	
	public static void printMK() {
		
		            String filePath = "logfiles/manual_output_" + Mutantigen.globalSimNum + ".txt";
        	String textToAppend = "About to save the .summary file.\n";

        try (FileWriter fw = new FileWriter(filePath, true)) {
            fw.write(textToAppend);
        } catch (IOException e) {
            System.err.println("An error occurred: " + e.getMessage());
        }
        
		try {
			String newFileName15 = "outfiles/out_" + Mutantigen.globalSimNum + ".summary";
			PrintStream summaryStream = new PrintStream(new FileOutputStream(newFileName15, true)); // append
			double sideBranchMut = (double) sideBranchMutations();
			double sideBranchOpp = sideBranchOpportunity();
			double sideBranchRate = sideBranchMut / sideBranchOpp;
			double trunkMut = (double) trunkMutations();
			double trunkOpp = trunkOpportunity();	
			double trunkRate = trunkMut / trunkOpp;		
			double mkRatio = trunkRate / sideBranchRate;
			summaryStream.printf("sideBranchRate\t%.4f\n", sideBranchRate);	
			summaryStream.printf("trunkRate\t%.4f\n", trunkRate);	
			summaryStream.printf("mkRatio\t%.4f\n", mkRatio);	
			summaryStream.close();
		} catch(IOException ex) {
			System.out.println("Could not write to file"); 
			System.exit(0);
		}
		
	}	
		
}
