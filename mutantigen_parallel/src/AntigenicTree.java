import java.util.*;
import cern.colt.matrix.*;
import java.io.*;

/**
 *
 * @author dar24
 */
public class AntigenicTree {
    
    	// fields
	private static Mutantigen root = null;	
	private static List<AntigenNode> nodes = new ArrayList<AntigenNode>();
        
        private static DoubleFactory2D factory2D = DoubleFactory2D.dense;
        private static DoubleMatrix2D distMatrix;
        
        private static ArrayList<ArrayList<Double>> distArray;
        
        private static ArrayList<Integer> treeTypes; 
        
        
        public static void initialize() {
            
            AntigenNode n = new AntigenNode(-1, 0.0);
            nodes.add(n);
            
            distArray = new ArrayList<ArrayList<Double>>();
            distArray.add(new ArrayList<Double>());
            distArray.get(0).add(0.0);
            
            //For new clean up methods
            treeTypes = new ArrayList<Integer>();
            treeTypes.add(0);
            
            if (Parameters.backgroundImmunity) {
                add(0, Parameters.backgroundDistance);
            }
            
            //distMatrix = factory2D.make(1,1);
            //distMatrix.setQuick(0,0,0.0);
        }
        
        public static int add(int parent, Double distance) {
            
            AntigenNode n = new AntigenNode(parent, distance);
            nodes.add(n);
            int count = treeTypes.size(); //nodes.size() - 1;
            updateDistArray(count, parent, distance);
            
            //For new clean up methods
            treeTypes.add(nodeCount()-1);
            
            //Old method with DoubleMatrix2D
            //updateDistMatrix(count);
            
            return (nodeCount() - 1); //count;
        }
        
        public static void remove(int type)  {
            
            int index = treeTypes.indexOf(type);
            if (index < 0) {
                System.out.println("WARNING: Cannot remove type from antigenic tree!");
            }
            treeTypes.remove(index);
            
            //Remove type from distance array
            distArray.remove(index);
            for (int j = 0; j < treeTypes.size(); j++) {
                distArray.get(j).remove(index);
            }
            
        }
        
        public static ArrayList<Integer> getTreeTypes()    {
            ArrayList<Integer> types = new ArrayList<Integer>();
            types.addAll(treeTypes);
            return types;
        }
        
        public static void updateDistArray(int count, int parent, double dist) {
            
            int parentIndex = treeTypes.indexOf(parent);
            
            // fill in new row
            distArray.add(new ArrayList<Double>());
            for(int i = 0; i < count; i++) { //for each column
                distArray.get(count).add(distArray.get(parentIndex).get(i) + dist);
            }
            distArray.get(count).add(0.0);
            
            // fill in new column
            for(int j = 0; j < count; j++) { //for each for
                distArray.get(j).add(distArray.get(j).get(parentIndex) + dist);
            }
            
        }
        
        public static ArrayList<ArrayList<Double>> getDistArray()
        {
            return distArray;
        }
        
        public static void updateDistMatrix(int count) {
            DoubleMatrix2D distCopy = (DoubleMatrix2D) distMatrix.copy();
            distMatrix = factory2D.make(count+1, count+1);
            for (int i = 0; i < count; i++) {
                for (int j = 0; j < count; j++) {
                    distMatrix.setQuick(i,j, distCopy.getQuick(i,j));
                }
            }
            // fill in perimeter
            double newDist;
            for (int i = 0; i < count+1; i++) {
                newDist = computeDistance(i,count);
                distMatrix.setQuick(i,count,newDist); // fill in last column
                distMatrix.setQuick(count,i,newDist); // fill in last row  
            }
            //System.out.println();
            
        }
        
        public static int nodeCount() {
            int count = nodes.size();
            return count;
        }
        
        public static double getDistance(int a1, int a2) {
            
            int a1Index = treeTypes.indexOf(a1);
            int a2Index = treeTypes.indexOf(a2);
            if (a1Index == -1 || a2Index == -1) {
                System.out.println("WARNING: Antigenic type not found in antigenic tree!");
            }
            
            double distance = distArray.get(a1Index).get(a2Index);
            
            //double distance = distMatrix.getQuick(a1,a2);
            
            return distance;
        }
        
        public static void printDistanceMatrix(ArrayList<Integer> types) {
            
            int nt = types.size();
            DoubleMatrix2D subMatrix = factory2D.make(nt,nt);
            subMatrix.setQuick(0,0,0.0);
            for (int i = 0; i < nt; i++) {
                for (int j = 0; j < nt; j++) {
                    int typeI = types.get(i);
                    int typeJ = types.get(j);
                    subMatrix.setQuick(i,j, computeDistance(typeI, typeJ));
                    //subMatrix.setQuick(i,j, distArray.get(typeI).get(typeJ));
                    //subMatrix.setQuick(i,j, distMatrix.getQuick(typeI, typeJ));
                }
            }
            
            try {
                String newFileName = "outfiles/out_" + Mutantigen.globalSimNum + ".antigenicDistances";
                File distFile = new File(newFileName);		
                distFile.delete();
                distFile.createNewFile();
                PrintStream stream = new PrintStream(distFile);
                
                for (int i = 0; i < nt; i++) {
                    for (int j = 0; j < nt; j++) {
                        stream.printf("%.4f", subMatrix.get(i,j));
                        if (j < nt - 1) {
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
        
        public static void printDistanceMatrixBetweenTips() {
            
            List<Virus> tips = VirusTree.getTips();
            int nt = tips.size();
            DoubleMatrix2D subMatrix = factory2D.make(nt,nt);
            subMatrix.setQuick(0,0,0.0);
            for (int i = 0; i < nt; i++) {
                for (int j = 0; j < nt; j++) {
                    Phenotype pI = tips.get(i).getPhenotype();
                    int typeI = pI.antigenicType();
                    Phenotype pJ = tips.get(j).getPhenotype();
                    int typeJ = pJ.antigenicType();
                    subMatrix.setQuick(i,j, computeDistance(typeI, typeJ));
                }
            }
            
            try {
                File distFile = new File("out.antigenicDistances");		
                distFile.delete();
                distFile.createNewFile();
                PrintStream stream = new PrintStream(distFile);
                
                for (int i = 0; i < nt; i++) {
                    for (int j = 0; j < nt; j++) {
                        stream.printf("%.4f", subMatrix.get(i,j));
                        if (j < nt - 1) {
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
        
        public static double computeDistance(int a1, int a2) {
            
            double distance = 0;
            int mrca = 0;
            if (a1 != a2) {
                ArrayList<Integer> ancestry1 = getAncestry(a1);
                ArrayList<Integer> ancestry2 = getAncestry(a2);
                boolean hitMRCA = false;
                int counter = 0;
                while (!hitMRCA) {
                    mrca = ancestry1.indexOf(ancestry2.get(counter));
                    if (mrca > -1) {
                        hitMRCA = true;
                    } else {
                        counter++;
                    }
                }
                double dist1 = sumDistances(ancestry1, mrca);
                double dist2 = sumDistances(ancestry2, counter);
                //if (ancestry1.size() > 2 & ancestry2.size() > 2) {
                    //System.out.println();
                //}
                distance = dist1 + dist2;
            }
            
            return distance;
        }
        
        public static ArrayList<Integer> getAncestry(int a) {
            ArrayList<Integer> ancestry = new ArrayList<Integer>();
            if (a > 0) {    
                ancestry.add(a);
                int parent = nodes.get(a).getParent();
                while (parent > 0) {
                    ancestry.add(parent);
                    parent = nodes.get(parent).getParent();
                }
            }
            ancestry.add(0);
            return ancestry;
        }
        
        public static double sumDistances(ArrayList<Integer> ancestry, int mrca) {
            double distance = 0;
            if (mrca > 0) {
                for (int i = 0; i < mrca; i++) {
                    distance += nodes.get(ancestry.get(i)).getDistance();
                }
            }
            return distance;
        }
        
        /**
        * New method for computing average host immunity
        * @param type antigenic type to compute immunity against
        * @param history antigenic types of previous infections
        * @return 
        */
        public static double getClosestDistance(int type, Phenotype[] history) {
            
            double closestDistance = 1.0; // was Double.POSITIVE_INFINITY before;
            if (history.length > 0) {
                for (int i = 0; i < history.length; i++) {
                    double thisDistance = getDistance(type, history[i].antigenicType());
                    if (thisDistance < closestDistance) {
                        closestDistance = thisDistance;
                    }
                }
            }
            return closestDistance;
            
        }
        
        
    
}
