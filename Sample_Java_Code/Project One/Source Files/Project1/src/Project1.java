import java.util.*;
import java.io.*;

/**
 *  The driver program for project1.  
 *  This driver reads in a file with matrixes with their order and then multiplies them.
 *
 *  @author  Jack Zhan
 *  @version 2016-02-20
 */


public class Project1 {
	
		//Storage for the Order of the Matrix
		public static List<Integer> LOrder = new ArrayList<Integer>();
		//Storage for the the Matrixes to be multiplied
		public static List<int[][]> LMatrix1 = new ArrayList<int[][]>();
		public static List<int[][]> LMatrix2 = new ArrayList<int[][]>();
		
		/**
	     *  Main entry point for the application.
	     */		
		public static void main (String args[]) {
			   
	        Project1          p        = new Project1();
	        MatrixMultiplier  mm       = new MatrixMultiplier();
	        String            fileName = "input.txt";
	       
	        System.out.println("Entered main() method.");
	        p.readInputFile(fileName);
	        mm.solve(p.LOrder,p.LMatrix1,p.LMatrix2);
	        return;
	    }
	
		/**
		 *  Opens, reads, and closes the file containing matrixes.
		*/

	private void readInputFile(String fileName) {
		
		int order = 0;
    	int[][] matrix1 = null;
    	int[][] matrix2 = null;
		String line = null;
		int flag  = 1;
		int flag2 = 0;
		int counter = 0;
		boolean Mswitch = true;
		
        System.out.println("Entered readInputFile() method.");
        
        try {
            // FileReader reads text files in the default encoding.
            FileReader fileReader = new FileReader(fileName);

            // Always wrap FileReader in BufferedReader.
            BufferedReader bufferedReader = new BufferedReader(fileReader);

            while((line = bufferedReader.readLine()) != null) {
                //Resets everything when there is a blank line and 
            	//adds the data to array lists
                if(line.isEmpty() || line.trim().equals("") || line.trim().equals("\n")){
                	flag  = 0;
                	flag2 = 0;
                	Mswitch = true;
                	LOrder.add(order);
                	LMatrix1.add(matrix1);
                	LMatrix2.add(matrix2);
                }
                if (flag==0) {
                	
                } else if(flag==1){
                	//Getting order and setting up the matrixes to be 
                	//added into the array list
                	order = Integer.parseInt(line);
                	matrix1 = new int[order][order];
                	matrix2 = new int[order][order];
                } else {
                	//Mswitch switches between Matrix1 and Matrix2
                	if (flag2 == order){
                		Mswitch = false;
                		flag2 = 0;
                	}
                	counter = 0;
                	if(Mswitch){
                		//Read in the data for the matrix
                		for(String temp : line.split("\\s")){
                			matrix1[flag2][counter] = Integer.parseInt(temp);
                			counter += 1;
                		}
                	} else {
                		for(String temp : line.split("\\s")){
                			matrix2[flag2][counter] = Integer.parseInt(temp);
                			counter += 1;
                		}
                	}
                	flag2 += 1;
                }
                flag += 1;
            }
        	LOrder.add(order);
        	LMatrix1.add(matrix1);
        	LMatrix2.add(matrix2);
 
            // Always close files.
            bufferedReader.close();         
        }
        catch(FileNotFoundException ex) {
            System.out.println(
                "Unable to open file '" + 
                fileName + "'");                
        }
        catch(IOException ex) {
            System.out.println(
                "Error reading file '" 
                + fileName + "'");                  
        }
        return;
    }

}
