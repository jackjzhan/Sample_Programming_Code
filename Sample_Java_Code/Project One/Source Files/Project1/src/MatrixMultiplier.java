import java.util.Arrays;
import java.util.List;

/**
 *  Multiplies two matrixes and calculates the total number of operations.
 * 
 *  @author  Jack Zhan
 *  @date    2016-02-20
 */

public class MatrixMultiplier {
    
	/**
     *  Create an object of the solver class.
     */

    public MatrixMultiplier() {
        System.out.println("Matrix Multiplier created.");
    }

	public double solve(List<Integer> LOrder, List<int[][]> LMatrix1, List<int[][]> LMatrix2) {
		
		int order = 0;
		int[][] matrix1 = null;
		int[][] matrix2 = null;
		int[][] matrix3 = null;
		int sum;
		int counter = 0;
		
		for( int index = 0; index<LOrder.size(); index++){
			// Grabbing the required inputs from Array list
			order = LOrder.get(index);
			matrix1 = LMatrix1.get(index);
			matrix2 = LMatrix2.get(index);
			matrix3 = new int[order][order];
			counter = 0;
			System.out.println("Multiplying Matrix of order " + order);
			//Algorithm for multiplying Matrixes
			for(int i=0; i<order; i++){
				for(int j=0; j<order; j++){
					sum=0;
					for(int k=0; k<order; k++){
						sum += matrix1[i][k]*matrix2[k][j];
						counter += 1;
					}
					matrix3[i][j] = sum;
				}
			}
			System.out.println("Matrix A = " + Arrays.deepToString(matrix1));
			System.out.println("Matrix B = " + Arrays.deepToString(matrix2));
			System.out.println("Solution = " + Arrays.deepToString(matrix3));
	        System.out.println("Total Number of operations = " + counter);

		}
        return 1.2;
	}
}
