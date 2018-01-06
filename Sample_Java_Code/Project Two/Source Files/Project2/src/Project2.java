import java.util.*;
import java.io.*;

/**
 *  The driver program for project2.  
 *  This driver reads in a file and hashes them into a table.
 *
 *  @author  Jack Zhan
 *  @version 2016-04-12
 */


public class Project2 {
	
		//Storage for the items
		public static List<Integer> Items = new ArrayList<Integer>();
		
		/**
	     *  Main entry point for the application.
	     */		
		public static void main (String args[]) {
			   
	        Project2          p        = new Project2();
	        HashTable         ht       = new HashTable();
	        HashTable2        ht2      = new HashTable2();
	        String            fileName = "input.txt";
	       
	        System.out.println("Entered main() method.");
	        p.readInputFile(fileName);
	        
	        ht.RunHashTable(Items, 120, 1, "Linear", 1, 120);
	        ht.RunHashTable(Items, 120, 1, "Quadratic", 1, 120);
	        ht2.RunHashTable(Items, 120, 1, 120);
	        
	        ht.RunHashTable(Items, 120, 1, "Linear", 1, 113);
	        ht.RunHashTable(Items, 120, 1, "Quadratic", 1, 113);
	        ht2.RunHashTable(Items, 120, 1, 113);
	        
	        ht.RunHashTable(Items, 40, 3, "Linear", 1, 41);
	        ht.RunHashTable(Items, 40, 3, "Quadratic", 1, 41);
	        
	        ht.RunHashTable(Items, 120, 1, "Linear", 2, 120);
	        ht.RunHashTable(Items, 120, 1, "Quadratic", 2, 120);
	        ht2.RunHashTable(Items, 120, 2, 120);
	        
	        return;
	    }
	
		/**
		 *  Opens, reads, and closes the file containing items.
		*/

	private void readInputFile(String fileName) {
		

		String line = null;
		
        System.out.println("Entered readInputFile() method.");
        
        try {
            // FileReader reads text files in the default encoding.
            FileReader fileReader = new FileReader(fileName);

            // Always wrap FileReader in BufferedReader.
            BufferedReader bufferedReader = new BufferedReader(fileReader);

            while((line = bufferedReader.readLine()) != null) {
                //Resets everything when there is a blank line and 
            	//adds the data to array lists
                if(line.isEmpty() ){
                } 
                else 
                {
                	//Getting order and setting up the matrixes to be 
                	//added into the array list
                	Items.add(Integer.parseInt(line));
                		
                } 
            }
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
