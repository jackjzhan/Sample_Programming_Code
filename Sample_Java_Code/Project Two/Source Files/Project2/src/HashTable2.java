import java.util.List;

/**
 *  Creates a Hash Table
 * 
 *  @author  Jack Zhan
 *  @date    2016-02-20
 */

public class HashTable2 {
    
	/**
     *  Create an object of the HashTable class.
     */
	
    private int MaxSize;       
    private HashEntry[] HTable;
	private int Function;
	private int Modulo;
	private int CollisionCount, CollisionNum;
    
	public HashTable2() 
    {
        CollisionCount = 0;
        CollisionNum = 0;
    }
    
    public void RunHashTable(List<Integer> Items, int size, int function, int modulo) 
    {
    	System.out.println("\nEntered Hash Table method.");
    	Function = function;
    	MaxSize = size;
    	Modulo = modulo;
    	System.out.println("Probe Type: Chaining Modulo: " + modulo);
    	HTable = new HashEntry[MaxSize];
        for (int i = 0; i < MaxSize; i++)
        {
        	HTable[i] = null;
        }
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
        if (HTable[HashVal] == null)
        {
        	HTable[HashVal] = new HashEntry(Key);
        	
        }
        else
        {
        	CollisionCount++;
        	CollisionNum++;
        	HashEntry Entry = HTable[HashVal];
        	while (Entry.getNext() != null)
        	{
        		CollisionNum++;
        		Entry = Entry.getNext();
        	}
        	Entry.setNext(new HashEntry(Key));
        }
    }
 
    /** Function to print HashTable **/
    private void printHashTable()
    {
    	int counter = 0;
    	boolean flag;
        System.out.println("\nHash Table: ");
        for (int i = 0; i < MaxSize; i++) 
        {
        	if (HTable[i] == null)
            {
        		if (counter == 4)
        		{
        			System.out.print("\n Index: " + i + " Value: Null");
        			counter = 0;
        		}
        		else
        		{
        			System.out.print(" Index: " + i + " Value: Null");
        			counter++;
        		}
            }
            else
            {
            	HashEntry Entry = HTable[i];
            	flag = false;
            	if (Entry.getNext() != null)
            	{
            		while (Entry.getNext() != null)
            		{
            			if (counter == 4)
            			{
            				System.out.print("\n Index: " + i + " Value: " + Entry.getKey());
            				counter = 0;
            				flag = true;
            			}
            			else
            			{
            				if (flag == false)
            				{
            					System.out.print(" Index: " + i + " Value: " + Entry.getKey());
            					counter++;
            					flag = true;
            				}
            				else
            				{
            					System.out.print(" Value: " + Entry.getKey());
            					counter++;
            				}
            			}
            			Entry = Entry.getNext();
            		}
            	}
            	else
            	{
            		if (counter == 4)
        			{
        				System.out.print("\n Index: " + i + " Value: " + Entry.getKey());
        				counter = 0;
        				flag = true;
        			}
        			else
        			{
        				System.out.print(" Index: " + i + " Value: " + Entry.getKey());
        				counter++;
        				flag = true;
        			}
            	}
            }
        }
        System.out.print("\nNumber of Collision: " + CollisionCount);
        System.out.print("\nTotal Number of Collision: " + CollisionNum);
    }
}