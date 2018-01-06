import java.util.*;
import java.io.*;

/**
 *  The driver program for project3.  
 *  This driver reads in a file and hashes them into a table.
 *
 *  @author  Jack Zhan
 *  @version 2016-05-03
 */


public class Project3 {
	
		//Storage for the items
		public static List<String> Items = new ArrayList<String>();
		private static String lcstring;
		
		/**
	     *  Main entry point for the application.
	     */		
		public static void main (String args[]) {
			   
	        Project3          p        = new Project3();
	        LCS               lcs      = new LCS();
	        String            fileName = "input.txt";
	       
	        System.out.println("Entered main() method.");
	        p.readInputFile(fileName);
	        for (String key1 : Items) 
	    	{
	        	for (String key2 : Items) 
		    	{
		        	if (key1 == key2)
		        	{
		        	}
		        	else
		        	{	
		        		System.out.println("LCS of " + key1 + " and " + key2 + " is:");
		        		lcstring = lcs.calLCS2(key1, key2);
		    	        System.out.println(lcstring);
		        	}
				}
			}
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
                if(line.isEmpty() )
                {
                } 
                else 
                {
                	//Getting order and setting up the matrixes to be 
                	//added into the array list
                	Items.add(line);
                		
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
