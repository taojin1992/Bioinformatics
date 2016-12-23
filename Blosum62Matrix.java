/*
***************************************************************************
*Jin Tao
*Student ID 11474660
Class Description:
Blosum62Matrix Class:
This class implements how to get the corresponding BLOSUM62 score 
for specific amino acid pairs from the BLOSUM62 scoring matrix. 
*****************************************************************************

* Blosum-62 substitution matrix
* #  A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V 
* A  4 -1 -2 -2  0 -1 -1  0 -2 -1 -1 -1 -1 -2 -1  1  0 -3 -2  0 
* R -1  5  0 -2 -3  1  0 -2  0 -3 -2  2 -1 -3 -2 -1 -1 -3 -2 -3 
* N -2  0  6  1 -3  0  0  0  1 -3 -3  0 -2 -3 -2  1  0 -4 -2 -3 
* D -2 -2  1  6 -3  0  2 -1 -1 -3 -4 -1 -3 -3 -1  0 -1 -4 -3 -3 
* C  0 -3 -3 -3  9 -3 -4 -3 -3 -1 -1 -3 -1 -2 -3 -1 -1 -2 -2 -1 
* Q -1  1  0  0 -3  5  2 -2  0 -3 -2  1  0 -3 -1  0 -1 -2 -1 -2 
* E -1  0  0  2 -4  2  5 -2  0 -3 -3  1 -2 -3 -1  0 -1 -3 -2 -2 
* G  0 -2  0 -1 -3 -2 -2  6 -2 -4 -4 -2 -3 -3 -2  0 -2 -2 -3 -3 
* H -2  0  1 -1 -3  0  0 -2  8 -3 -3 -1 -2 -1 -2 -1 -2 -2  2 -3 
* I -1 -3 -3 -3 -1 -3 -3 -4 -3  4  2 -3  1  0 -3 -2 -1 -3 -1  3 
* L -1 -2 -3 -4 -1 -2 -3 -4 -3  2  4 -2  2  0 -3 -2 -1 -2 -1  1 
* K -1  2  0 -1 -3  1  1 -2 -1 -3 -2  5 -1 -3 -1  0 -1 -3 -2 -2 
* M -1 -1 -2 -3 -1  0 -2 -3 -2  1  2 -1  5  0 -2 -1 -1 -1 -1  1 
* F -2 -3 -3 -3 -2 -3 -3 -3 -1  0  0 -3  0  6 -4 -2 -2  1  3 -1 
* P -1 -2 -2 -1 -3 -1 -1 -2 -2 -3 -3 -1 -2 -4  7 -1 -1 -4 -3 -2 
* S  1 -1  1  0 -1  0  0  0 -1 -2 -2  0 -1 -2 -1  4  1 -3 -2 -2 
* T  0 -1  0 -1 -1 -1 -1 -2 -2 -1 -1 -1 -1 -2 -1  1  5 -2 -2  0 
* W -3 -3 -4 -4 -2 -2 -3 -2 -2 -3 -2 -3 -1  1 -4 -3 -2 11  2 -3 
* Y -2 -2 -2 -3 -2 -1 -2 -3  2 -1 -1 -2 -1  3 -3 -2 -2  2  7 -1 
* V  0 -3 -3 -3 -1 -2 -2 -3 -3  3  1 -2  1 -1 -2 -2  0 -3 -1  4 
*/
public class Blosum62Matrix {
	//Initially, I define a constant class variable matrix, which is the BLOSUM62 scoring matrix
	private static final int[][] matrix = 
		{
			{ 4, -1, -2, -2,  0, -1, -1,  0, -2, -1, -1, -1, -1, -2, -1,  1,  0, -3, -2,  0},
			{-1,  5,  0, -2, -3,  1,  0, -2,  0, -3, -2,  2, -1, -3, -2, -1, -1, -3, -2, -3},
			{-2,  0,  6,  1, -3,  0,  0,  0,  1, -3, -3,  0, -2, -3, -2,  1,  0, -4, -2, -3},
			{-2, -2,  1,  6, -3,  0,  2, -1, -1, -3, -4, -1, -3, -3, -1,  0, -1, -4, -3, -3},
			{ 0, -3, -3, -3,  9, -3, -4, -3, -3, -1, -1, -3, -1, -2, -3, -1, -1, -2, -2, -1},
			{-1,  1,  0,  0, -3,  5,  2, -2,  0, -3, -2,  1,  0, -3, -1,  0, -1, -2, -1, -2},
			{-1,  0,  0,  2, -4,  2,  5, -2,  0, -3, -3,  1, -2, -3, -1,  0, -1, -3, -2, -2},
			{ 0, -2,  0, -1, -3, -2, -2,  6, -2, -4, -4, -2, -3, -3, -2,  0, -2, -2, -3, -3},
			{-2,  0,  1, -1, -3,  0,  0, -2,  8, -3, -3, -1, -2, -1, -2, -1, -2, -2,  2, -3},
			{-1, -3, -3, -3, -1, -3, -3, -4, -3,  4,  2, -3,  1,  0, -3, -2, -1, -3, -1,  3},
			{-1, -2, -3, -4, -1, -2, -3, -4, -3,  2,  4, -2,  2,  0, -3, -2, -1, -2, -1,  1},
			{-1,  2,  0, -1, -3,  1,  1, -2, -1, -3, -2,  5, -1, -3, -1,  0, -1, -3, -2, -2},
			{-1, -1, -2, -3, -1,  0, -2, -3, -2,  1,  2, -1,  5,  0, -2, -1, -1, -1, -1,  1},
			{-2, -3, -3, -3, -2, -3, -3, -3, -1,  0,  0, -3,  0,  6, -4, -2, -2,  1,  3, -1},
			{-1, -2, -2, -1, -3, -1, -1, -2, -2, -3, -3, -1, -2, -4,  7, -1, -1, -4, -3, -2},
			{ 1, -1,  1,  0, -1,  0,  0,  0, -1, -2, -2,  0, -1, -2, -1,  4,  1, -3, -2, -2},
			{ 0, -1,  0, -1, -1, -1, -1, -2, -2, -1, -1, -1, -1, -2, -1,  1,  5, -2, -2,  0},
			{-3, -3, -4, -4, -2, -2, -3, -2, -2, -3, -2, -3, -1,  1, -4, -3, -2, 11,  2, -3},
			{-2, -2, -2, -3, -2, -1, -2, -3,  2, -1, -1, -2, -1,  3, -3, -2, -2,  2,  7, -1},
			{ 0, -3, -3, -3, -1, -2, -2, -3, -3,  3,  1, -2,  1, -1, -2, -2,  0, -3, -1,  4}
		};
	/*private static int getIndex(char input): This method returns the column or row index  
	for input character representing specific amino acid in the BLOSUM62 scoring matrix.*/
	 private static int getIndex(char input) {
		 int out=0;
		 switch((String.valueOf(input)).toUpperCase().charAt(0)) {
		 	case 'A': 
		 		out=0;
		 		break;
	    	case 'R': 
	    		out=1;
		 		break;
	    	case 'N': 
	    		out=2;
		 		break;
	    	case 'D': 
	    		out=3;
		 		break;
	    	case 'C': 
	    		out=4;
		 		break;
	    	case 'Q': 
	    		out=5;
		 		break;
	    	case 'E': 
	    		out=6;
		 		break;
	    	case 'G': 
	    		out=7;
		 		break;
	    	case 'H': 
	    		out=8;
		 		break;
	    	case 'I': 
	    		out=9;
		 		break;
	    	case 'L': 
	    		out=10;
		 		break;
	    	case 'K': 
	    		out=11;
		 		break;
	    	case 'M': 
	    		out=12;
		 		break;
	    	case 'F': 
	    		out=13;
		 		break;
	    	case 'P':
	    		out=14;
		 		break;
	    	case 'S': 
	    		out=15;
		 		break;
	    	case 'T': 
	    		out=16;
		 		break;
	    	case 'W': 
	    		out=17;
		 		break;
	    	case 'Y': 
	    		out=18;
		 		break;
	    	case 'V': 
	    		out=19;
		 		break;
		 }
		return out;
	 }
	 /*private static int getScore(int row, int col): This method returns the corresponding BLOSUM62 score 
	  * given the position in the BLOSUM62 scoring matrix.
	  */
	 private static int getScore(int row, int col) {
			 return matrix[row][col];
	 }
	 /*public static int getDiagScore(char a1, char a2): This methods returns the corresponding BLOSUM62 score 
	  * for specific amino acid pairs from the BLOSUM62 scoring matrix. As a class method, 
	  * in the global sequence alignment we directly use Blosum62Matrix.getDiagScore(char a1, char a2) 
	  * to get the corresponding BLOSUM62 score for specific amino acid pairs.
	  */
	 public static int getDiagScore(char a1, char a2) {
	    	return getScore(getIndex(a1), getIndex(a2));  	
	    }
}
