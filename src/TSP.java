import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.Arrays;
import java.util.BitSet;
import java.util.HashMap;
import java.util.Map;
import java.util.StringTokenizer;

/** Class: TSP.java
 *  @author Yury Park
 *
 *  This class - Algorithm for solving the TSP problem in O(n^2 * 2^n) time (as opposed to O(n!) via naive brute-force search)
 *  
 *  Purpose - 
 *  
 *  You are given a data file, whose first line indicates the number of cities. 
 *  Each city is a point in the plane, and each subsequent line indicates the x- and y-coordinates of a single city.
 *  
 *  The distance between two cities is defined as the Euclidean distance --- that is, two cities at locations 
 *  (x,y) and (z,w) have distance sqrt[(x-z)^2 + (y-w)^2] between them.
 *  
 *  This class figures out the minimum cost of a traveling salesman tour for this instance starting from city 0.
 *  
 *  =======================================================================================
 *  NOTE: This algorithm my take up a lot of memory depending on the size of the data file.
 *  The following VM argument may help:
 *  -Xmx4096m -Xms4096m -Xmn256m -Xss16m
 *  =======================================================================================
 */
public class TSP {
	private String filename;
	private Map<Integer, Vertex> vertexMap;	//maps an int value (label no.) to the associated city vertex.
	
	/* The below array will be implemented in the following format: A[S][j] = some floating number, where:
	 * 
	 * S: given subset of Vertices
	 * j: given Vertex.
	 * 
	 * This array keeps track of the MIN. distance from starting Vertex (assumed to be Vertex 0) to some
	 * destination vertex j (where j may or may not be the same as the starting Vertex), PROVIDED that
	 * the solution path MUST visit each vertex in the given subset of vertices S EXACTLY once.
	 * */
	private float[][] A;			
	private float[][] distances;	//Caches the distance between every pair of Vertices to avoid re-computing them
	private int N;					//The total no. of vertices
	private static boolean debugOn;
	
	/**
	 * 2-arg constructor.
	 * @param filename
	 * @param debug
	 */
	public TSP (String filename, boolean debug) {
		this.filename = filename;
		this.vertexMap = new HashMap<Integer, Vertex>();	//This will map city vertex label no. to the actual city Vertex object.
		debugOn = debug;
		build();	//custom method invoke
	}
	
	 /**
     * Method: build
     * Reads in the given data text file and constructs a Graph object composed of 
     * a bunch of Vertices (representing cities on a map).
     */
    private void build(){  
        try {
            BufferedReader rd = new BufferedReader(new FileReader(new File(filename)));
            
            //Read the first line, which contains info re: total no. of vertices
            String line = rd.readLine();
            StringTokenizer tokenizer = new StringTokenizer(line);  
            this.N = Integer.parseInt(tokenizer.nextToken());	//N = total no. of vertices
            this.distances = new float[N][N];					//Initialize array to keep track of distances between every pair of vertices
            
            /* Initialize the array with size of [2^(N-1)][N]. We need 2^(N-1) for the first dimension
             * because we'll be saving the integer value of EVERY possible binary permutation consisting of N-1 digits.
             * For example, if there are N = 4 vertices overall, then we need every binary permutation consisting of
             * 4-1 = 3 digits, as follows:
             * 
             * 000 001 010 100 011 101 110 111
             * 
             * The above is a total of 2^(N-1) = 2^(4-1) = 2^3 = 8 permutations. What do the permutations represent?
             * Each permutation represents a possible subset of vertices that should be explored (where 0 = unexplored, 1 = explored)
             * at any given time. Read from right to left, each digit represents Vertex index 1 thru Vertex index N-1.
             * (We designate Vertex index 0 as the starting Vertex, which is always explored so we omit it from the
             * above permutation to save memory.)
             * 
             * We need N for the 2nd dimension of the array because there is a total of N vertices, indexed from 0 to N-1.
             * */
            this.A = new float[(int)Math.pow(2, N-1)][N];
            
            //Pre-fill every slot in the array with Infinity.
            for (int i = 0; i < A.length; i++) {
            	float[] tempArr = A[i];
            	Arrays.fill(tempArr, Float.MAX_VALUE);	//Arrays.fill method only work on 1-dimensional arrays, thus the for loop
            }
            
            for (int j = 0; j < N; j++) {	//read in all n vertices (let's label them from 0 to N-1)
            	line = rd.readLine();
                tokenizer = new StringTokenizer(line);
                
                /* Each line contains the (x, y) location of each succeeding vertex. */
                float x = Float.parseFloat(tokenizer.nextToken());
                float y = Float.parseFloat(tokenizer.nextToken());
                Vertex v = new Vertex(j, x, y);	//Construct vertex with the given label no. and (x, y) position.
                
				this.vertexMap.put(j, v);	//Add vertex to the HashMap. the calculateDistance() method below uses this Map.
				
				/* Pre-compute and cache the distances between every pair of vertices. Do this as we're reading in data from
				 * txt file. 
				 * 
				 * Note that there is an implied "Edge" connecting EVERY possible pair of vertices. So this is a Clique graph,
				 * which by definition is a graph where EVERY vertex neighbors EVERY other vertex. */
				for (int i = 0; i < j; i++) {
					this.distances[i][j] = calculateDistance(i, j);
					this.distances[j][i] = this.distances[i][j];	//Distance from j to i = distance from i to j.
				}
            }  
            //end for i
            
            if (debugOn) {
            	System.out.println("Finished reading from file. Distance Adjacency Matrix:\n{");
            	for (int i = 0; i < N; i++) {
            		for (int j = 0; j < N; j++) {
//            			System.out.print((j == 0 ? "{" : "") + this.distances[i][j] + (j == N-1 ? "}," : ", "));
            			System.out.printf("%s%9f%s", j == 0 ? " {" : "", this.distances[i][j], j == N-1 ? i == N-1 ? "}" : "}," : ", ");
            		}
            		System.out.println();
            	}
            	System.out.println("}");
            }
            rd.close();	//close the reader when done.
            this.vertexMap = null;	//We no longer need this HashMap since we'll no longer invoke calculateDistance() method. 
            
        } catch (FileNotFoundException e) {   
            e.printStackTrace();  
        } catch (IOException e) {  
            e.printStackTrace();  
        }  
    }  
    //end private void build()
    
    /**
     * Method: calculateDistance
     * @param v1 label no. of a vertex
     * @param v2 label no. of a vertex
     * @return the Euclidean distance between v1 and v2.
     */
    private float calculateDistance(int v1, int v2) {
    	Vertex vertex1 = this.vertexMap.get(v1);
    	Vertex vertex2 = this.vertexMap.get(v2);
    	
    	/* Use distance formula. 
    	 * Two cities at locations (x1,y1) and (x2,y2) have distance sqrt[(x1-x2)^2 + (y1-y2)^2] between them */
    	return (float)Math.sqrt((vertex1.x - vertex2.x) * (vertex1.x - vertex2.x) +
    							(vertex1.y - vertex2.y) * (vertex1.y - vertex2.y));
    }
    
    /**
     * Method: compute
     *         Computes the total distance of the most efficient path for the following task:
     *         start from Vertex 0 and visit every vertex EXACTLY once, then return to the starting vertex.
     *         (NOTE: For convenience we designated the starting vertex to be vertex 0, but the optimal
     *          distance will stay the same no matter which starting vertex we pick.)
     * 
     * @return the min. distance needed to visit every vertex EXACTLY once and to return to the starting vertex. 
     */
    private float compute() {
    	
    	/* BIG PICTURE: We will employ roughly the following optimal substructure:
    	 * 
    	 * For every destination city vertex, j, that is an element of the set {0, 1, 2, ... n-1}, and 
    	 * for every subset S of the aforementioned set that contains both 0 and j, 
    	 * we let A[S][j] = minimum length of a path from 0 to j that visits only the vertices contained in S
    	 * EXACTLY ONCE each. 
    	 * 
    	 * ===================================================================================================
    	 * 
    	 * HOW WE COMPUTE THE RECURRENCE:
    	 * 
    	 * Let P be a shortest path from 0 to j that visits all the vertices in S exactly once each.
    	 * (here we'll assume that |S| >= 2.) Then, if the last "hop" in the path P is k -> j, where k
    	 * is some vertex in S, then we know that P' is a shortest path from 0 to k that visits all the
    	 * vertices in S - {j} exactly once each.
    	 * 
    	 * CORRESPONDING RECURRENCE:
    	 * A[S][j] = min_{for all k in S, where k != j} (A[S-{j}][k] + distance[k][j])
    	 * 
    	 * ===================================================================================================
    	 * PSEUDOCODE FOR DYN. PROGRAMMING ALGORITHM:
    	 * 
    	 * A[S][0] = 0 if S = {0}, else infinity
    	 * 
    	 * For m = 1 to N-1:  //where m = size of subproblem
    	 *     For each subset S of the given set {0, 1, ..., n-1} of cardinality m that contains vertex 0:
    	 *         For each element j in S where j != starting vertex 0:
    	 *             A[S][j] = (See the CORRESPONDING RECURRENCE formula above)
    	 * 
    	 * //Return the min. cost of travelling from vertex 0 to each possible value of j (1 thru N-1) 
    	 * //while visiting every vertex once, PLUS
    	 * //the distance of the final hop from j back to the starting vertex 0.
    	 * Return min_{j=1 to N-1} ( A[{0, 1, ..., N-1}][j] + distance[j][0] )
    	 * */
    	
    	/* Base case. Min. distance between starting vertex and vertex 0 (AKA from vertex 0 to itself) is zero
    	 * IFF the set of vertices that must be explored EXACTLY once in the solution path consists of vertex 0 only.
    	 * (otherwise, the distance is infinity. We don't have to assign infinity in this method because we 
    	 * already did that during the initial text data file reading.) 
    	 **/
    	A[0][0] = 0f;	
    	
    	//Now go thru every other vertex (from Vertex index 1 to N-1) 
    	for (int m = 1; m < N; m++) {
    		System.out.printf("========================================\nNow starting iteration m = %s...\n", m);
    		
    		/* Gosper's hack. (Look it up on the web) 
    		 * This creates every possible permutation of binary numbers where m = total number of 1's
    		 * and N = total number of digits. So for example, if N = 5 and m = 3, then we'd have the following permutation:
    		 * 
    		 * 00111 (integer value: 7)
    		 * 01011 (integer value: 11)
    		 * 01101 (integer value: 13)
    		 * 01110 (integer value: 14)
    		 * 10011 ...and so on...
    		 * 10101 
    		 * 10110 
    		 * 11001 
    		 * 11010 
    		 * 11100 
    		 * 
    		 * See my comments in the build() method above to see what these permutations represent. 
    		 * 
    		 * NOTE: Slight modification to save memory...Since the startingVertex is 0 and that is ALWAYS visited,
    		 * it is redundant to include vertex 0 as the rightmost bit in the above permutation (because then, the rightmost
    		 * bit would ALWAYS be 1). So we simply omit Vertex 0 from the above permutations and assume that it's
    		 * always visited. */
    		for (int s = (1 << m) - 1; (s >>> N-1) == 0; s = nextCombo(s)) {	//nextCombo(s) invokes Gosper's hack method.
    			int setOfAllVerticesToBeExplored = s;
    			
    			//convert to binary
    			String setOfAllVerticesToBeExplored_Binary = Integer.toBinaryString(setOfAllVerticesToBeExplored);
    			if (debugOn) 
    				System.out.printf("setOfAllVerticesToBeExplored:\t%s (binary: %s) <-- This is the set of vertices that must be explored"
    				       	        + " EXACTLY once (starting vertex 0 is omitted to save memory, but it's ALWAYS assumed to be explored)\n"
    				       	        + "\t\t\t\t\t\tHow to read the binary expression: Read it from right to left. Assign the rightmost digit to vertex 1,\n"
    				       	        + "\t\t\t\t\t\tthe 2nd rightmost digit to vertex 2, and so on. If a digit equals 1, that vertex must be explored EXACTLY once.\n"
    				       	        + "\t\t\t\t\t\tIf a digit equals 0, that vertex must not be explored at this time.\n",
    				       	        setOfAllVerticesToBeExplored, Integer.toBinaryString(setOfAllVerticesToBeExplored));
    			
    			int vertexOtherThanStartingVertex = 1;	//Since starting Vertex label is 0, we begin with 1.
    			
    			/* Now we traverse the binary string from right to left. We basically want to see which vertices 
    			 * in this binary string are set as 1, AKA "should be explored EXACTLY ONCE", and then for every 
    			 * such vertex, we will invoke the getMin() method to get the min. distance. */
    			for (int j = setOfAllVerticesToBeExplored_Binary.length() - 1; j >= 0; j--, vertexOtherThanStartingVertex++) {
    				/* If a given position is 0, then it means the vertex at this position is NOT explored. So move on. */
    				if (setOfAllVerticesToBeExplored_Binary.charAt(j) == '0') continue;
    				A[setOfAllVerticesToBeExplored][vertexOtherThanStartingVertex] =
    						getMin(setOfAllVerticesToBeExplored, vertexOtherThanStartingVertex);	//custom method!

    				if (debugOn)
    					System.out.printf("A[%s][%s] (Min. distance from 0 to %s that includes the subset of vertices that must be explored EXACTLY once"
    							+ ", where %s = %s in binary) = %s\n",
    							setOfAllVerticesToBeExplored,
    							vertexOtherThanStartingVertex,
    							vertexOtherThanStartingVertex,
    							setOfAllVerticesToBeExplored,
    							setOfAllVerticesToBeExplored_Binary,
    							A[setOfAllVerticesToBeExplored][vertexOtherThanStartingVertex]);
    			}
    			//end for j
    		}
    		//end for s
    	}
    	//end for m
    		
    	if (debugOn) {
    		System.out.println("\nNow computing the most efficient path length to traverse the entire graph, starting at vertex 0...");
    	}
    	
    	float ret = Float.MAX_VALUE;
    	BitSet bitset = new BitSet();	//BitSet is a built-in java class.
    	
    	/* For this final step, we want a binary string where every vertex is visited exactly once.
    	 * For example, if N = 5, the string is: 
    	 * 
    	 * 1111 
    	 * 
    	 * where we exclude the starting vertex from the above string (to save array memory, and because
    	 * the starting vertex is ALWAYS explored anyway, so it's redundant to include it).
    	 * 
    	 * So, we would set every index in the above binary string (from index 0 to 3) to 1.
    	 * Note that for an arbitrary N, this would be basically setting (from index 0 to N-2) to 1.
    	 * 
    	 * Afterwards, we want the integer value of the binary string, which in the above example case is:
    	 * 1 + 2 + 4 + 8 = 15. */
    	for (int i = 0; i <= this.N - 2; i++) bitset.set(i, true);	
    	int bitsetInt = this.convertToInt(bitset);	//custom method to convert binary String to integer
    	
    	/* Now we go from Vertex index label 1 thru N-1. */
    	for (int j = 1; j < this.N; j++) {
    		float distance = A[bitsetInt][j];
    		if (distance == Float.MAX_VALUE) continue;
    		ret = Math.min(ret, distance + this.distances[j][0]);
    	}
    	
    	return ret;
    }
    //end private float compute
    
    /**
     * Method: nextCombo
     *         This is Gosper's Hack. Look it up on the web. Also see my comments in the compute() method above.
     * @param x given integer.
     * @return the next larger permutation of the binary string, expressed in integer form.
     */
    private int nextCombo (int x) {
    	// moves to the next combination with the same number of 1 bits
    	int u = x & (-x);
    	int v = u + x;
    	return v + (((v ^ x) / u) >> 2);
    }
	
    /**
     * Method: getMin
     * @param s Integer representation of a subset of vertices that must be explored EXACTLY once.
     * @param j index no. of a vertex.
     * @return Goes thru every neighbor of Vertex j, then computes A[sMinusj][k] + distance between k and j,
     *         and returns the minimum of those values.
     */
    private float getMin(int s, int j) {
    	float ret = Float.MAX_VALUE;
    	int sMinusj = s - (int)Math.pow(2, j-1);	//This computes S - {j}.
    	
    	//Get min(A[sMinusj][k] + Ckj) (where k is j's neighbor)
    	for (int k = 0; k < this.N; k++) {
    		if (k == j) continue;	//if k == j, then it isn't j's neighbor, so moving on. In every other case, k is j's neighbor.
    		float distance =  A[sMinusj][k];
    		if (distance == Float.MAX_VALUE) continue;
    		ret = Math.min(ret, distance + this.distances[k][j]);
    	}
    	return ret;
    }
   
    /**
     * Method: convertToInt
     * @param bs a BitSet object. (representing a binary string)
     * @return the Integer value of the given BitSet.
     */
 	private int convertToInt(BitSet bs) {
		int ret = 0;
		for (int i = 0; i < bs.length(); i++) {
			if (bs.get(i) == true) {
				ret += (int)Math.pow(2, i);
			}
		}
		return ret;
	} 
    
    /**
	 * Class: Vertex
	 *        Vertex object. Inner (nested) class. 
	 */
	private class Vertex {
		int lbl;	//name (label no.) of vertex
		float x, y;	//the (x, y) location of this Vertex on a 2-D plane.
		
		
		/**
		 * 1-arg constructor.
		 * @param lbl the name (label) to assign to this vertex
		 */
		Vertex(int lbl, float x, float y) {
			this.lbl = lbl;
			this.x = x;
			this.y = y;
		}
		
		@Override
		/**
		 * Method: toString
		 * @return the label no. of this Vertex.
		 */
		public String toString() {
			return String.format("%s (%s, %s)", lbl, x, y);
		}
	}
	//end private class Vertex

	/**
	 * Method: printResults
	 * @param file
	 * @param debugOn
	 */
	public static void printResults(String file, boolean debugOn) {
		long startTime = System.currentTimeMillis();
		System.out.printf("============================================\nNow running file %s...\n", file);
		TSP tsp = new TSP(file, debugOn);
		System.out.printf("--------------------------------------------------------------------------------------------------\n"
				          + "Min. distance to visit every vertex EXACTLY once and to return to the starting vertex: %s\n", tsp.compute());
		System.out.printf("Time elapsed (in millisecs): %s\n", System.currentTimeMillis() - startTime);
	}
	
	/**
	 * Method: main
	 * @param args
	 */
	public static void main(String... args) { 
		if (args.length == 0) {
			System.out.println("Please input a filename. Ending program...");
		}
		
		/* End user will input something like "tsp_BIG.txt" (without the quotes). */
		for (String fileName : args) {			
			printResults(fileName, false);  //run method with debug set to true or false
		}
	}
	//end main
}