import java.util.List;

/**
 *  Creates a Hash Table
 * 
 *  @author  Jack Zhan
 *  @date    2016-02-20
 */

public class HashTable {
    
	/**
     *  Create an object of the HashTable class.
     */
	
    private int MaxSize;       
    private int[][] HTable;
    private String Probe;
	private int Function;
	private int Bucket;
	private int Modulo;
	private int CollisionCount, CollisionNum, InvalidData;
    private int c1 = 3;
    private int c2 = 1;
	public HashTable() 
    {
        CollisionCount = 0;
        CollisionNum   = 0;
        InvalidData    = 0;
    }
    
    public void RunHashTable(List<Integer> Items, int size, int bucket, String probe, int function, int modulo) 
    {
    	System.out.println("\nEntered Hash Table method.");
    	Probe    =  probe;
    	Function = function;
    	Bucket   = bucket;
    	MaxSize  = size;
    	Modulo   = modulo;
    	System.out.println("Bucket Size: " + Bucket + " Probe Type: " + Probe + " Modulo: " + modulo);
        HTable = new int[MaxSize][Bucket];
    	for (int key : Items) 
    	{
			insert(key, hash(key));
		}
        printHashTable();
    }
    
    /** Function to get hash code of a given Key **/
    private int hash(int Key) 
    {
        if (Function == 1)
        {
        	return Key % Modulo;
        } 
        else if (Function == 2)
        {
        	return (Key^2) % Modulo;
        }
        else
        {
        	System.out.println("Invalid Value for Hash Function.");
        	return 0;
        }
    }    
 
    /** Function to insert Key-value pair **/
    private void insert(int Key, int Hash) 
    {                
        int HashVal = Hash;
        int value;
        boolean flag = false;
        for(int i=0; i<MaxSize; i++)
        {
            if (Probe == "Linear")
            {
            	value = (HashVal + i) % MaxSize;
            }
            else if (Probe == "Quadratic")
            {
            	value = (HashVal + c1*i + c2*i*i) % MaxSize;
            }
            else
            {
            	System.out.println("Invalid Value for Probe.");
            	return;
            }
            for (int k=0; k<Bucket; k++)
            {
            	if (HTable[value][k] == 0)
            	{
            		HTable[value][k] = Key;
            		return;
            	}
            	else
            	{
                	if (flag == false)
                	{
                		flag = true;
                		CollisionCount++;
                	}
                	CollisionNum++;
            	}
            }
        }
        InvalidData++;
        return;
    }
 
    /** Function to print HashTable **/
    private void printHashTable()
    {
    	int counter = 0;
    	boolean flag;
    	int count;
        System.out.println("\nHash Table: ");
        if (Bucket == 1)
		{
        	count=4;
		}
        else
        {
        	count=5;
        }
        
        for (int i = 0; i < MaxSize; i++) 
        {
        	flag = false;
        	for (int j = 0; j < Bucket; j++)
			{
        		if (counter == count)
        		{
        			if (flag==false)
        			{
        				System.out.println(" Index: " + i + " Value: " + HTable[i][j]);
        				flag = true;
        			}
        			else
        			{
        				System.out.println(" Value: " + HTable[i][j]);
        			}
        		}
        		else
        		{
        			if (flag==false)
        			{
        				System.out.print(" Index: " + i + " Value: " + HTable[i][j]);
        				flag = true;
        			}
        			else
            		{
            			System.out.print(" Value: " + HTable[i][j]);
            		}
        		}
        		counter++;
        		if (counter == count+1)
        		{
        			counter = 0;
        		}
			}
        }
        System.out.println("Number of Collision: " + CollisionCount);
        System.out.println("Total Number of Collision: " + CollisionNum);
        System.out.println("Data not inputed into Hash Table: " + InvalidData);
    }
}