/*
**************************************************************************************************************************
*Jin Tao
*Student ID 11474660
Class Description:
GlobalSequenceAlignmentTable Class: 
This class implements the Global Sequence Alignment Dynamic Programming Algorithm for DNA sequence and protein sequence. 
The essential frameworks for two kinds of sequence alignment are the same except that we follow the BLOSUM62 scoring matrix 
to calculate the score of matching one amino acid with another while we use the constant score of matching/mismatching, 
namely, 1 and -1. 
***************************************************************************************************************************
*/
import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.Scanner;

public class GlobalSequenceAlignmentTable { 
	/*
	 * protected static void initializeScores(Cell[][] scoreTable): 
	 * This method initializes the scores in all cells in the first row and the first column in the alignment matrix. 
	 */
	protected static void initializeScores(Cell[][] scoreTable) {	
		int gapPenalty=2;
		scoreTable[0][0].setScore(0);
		for(int i=1;i<scoreTable[1].length;i++) {
			scoreTable[0][i].setScore(scoreTable[0][i-1].getScore()-gapPenalty);
		}
		for(int i=1;i<scoreTable.length;i++) {
			scoreTable[i][0].setScore(scoreTable[i-1][0].getScore()-gapPenalty);
		}
	}
	/*
	 * protected static void initializePointers(Cell[][] scoreTable): 
	 * This method initializes the arrows in all cells in the first row and the first column in the alignment matrix. 
	 */
	protected static void initializePointers(Cell[][] scoreTable) {
		scoreTable[0][0].setPrevCell(null);
		
		for(int i=1;i<scoreTable[0].length;i++) {
			scoreTable[0][i].setPrevCell(scoreTable[0][i-1]);
		}
		for(int i=1;i<scoreTable.length;i++) {
			scoreTable[i][0].setPrevCell(scoreTable[i-1][0]);
		}
	}
	/*
	 * protected static void initialize(Cell[][] scoreTable): 
	 * This method builds each cell inside the alignment matrix and initialize the whole alignment matrix. 
	 */
	protected static void initialize(Cell[][] scoreTable) {
		for(int i=0;i<scoreTable.length;i++) {//row
			for(int j=0;j<scoreTable[i].length;j++) {
				scoreTable[i][j]=new Cell(i,j);
			}
		}
		initializeScores(scoreTable);
		initializePointers(scoreTable);
	}
	/*
	 * protected static void fillInCellDNASeq(String sequence1,String sequence2, Cell currentCell,Cell cellAbove, Cell cellToLeft, Cell cellAboveLeft):
	 *  This method fills in one Cell object in the alignment matrix for the DNA sequence alignment. 
	 *  For DNA sequence alignment, we use 1 and -1 as the score of matching/mismatching; use -2 as the gap penalty.
	 */
	protected static void fillInCellDNASeq(String sequence1,String sequence2, Cell currentCell,Cell cellAbove, Cell cellToLeft, Cell cellAboveLeft) {
		int gapPenalty=2;
		int matchScore=1;
		int verticalScore=cellAbove.getScore()-gapPenalty;
		int horizontalScore=cellToLeft.getScore()-gapPenalty;
		int diagonalScore=cellAboveLeft.getScore();
		if(sequence2.charAt(currentCell.getRow()-1)==sequence1.charAt(currentCell.getCol()-1)) {
			diagonalScore=diagonalScore+matchScore;
		}
		else {
			diagonalScore=diagonalScore-matchScore;
		}
		if(verticalScore>=horizontalScore) {
			if(verticalScore>=diagonalScore) {
				currentCell.setScore(verticalScore);
				currentCell.setPrevCell(cellAbove);
			}
			else {
				currentCell.setScore(diagonalScore);
				currentCell.setPrevCell(cellAboveLeft);
			}
		}
		else {
			if(horizontalScore>=diagonalScore) {
				currentCell.setScore(horizontalScore);
				currentCell.setPrevCell(cellToLeft);
			}
			else {
				currentCell.setScore(diagonalScore);
				currentCell.setPrevCell(cellAboveLeft);
			}
		}
	}
	/*
	 * protected static void fillInCellProSeq(String sequence1,String sequence2, Cell currentCell,Cell cellAbove, Cell cellToLeft, Cell cellAboveLeft): 
	 * This method fills in one Cell object in the alignment matrix for the protein sequence alignment. 
	 * For protein sequence alignment, we follow the BLOSUM62 scoring matrix to calculate the score of matching one amino acid with another; 
	 * use -2 as the gap penalty.
	 */	
	protected static void fillInCellProSeq(String sequence1,String sequence2, Cell currentCell,Cell cellAbove, Cell cellToLeft, Cell cellAboveLeft) {
		int gapPenalty=2;
		int matchScore=Blosum62Matrix.getDiagScore(sequence1.charAt(currentCell.getCol()-1), sequence2.charAt(currentCell.getRow()-1));
		int verticalScore=cellAbove.getScore()-gapPenalty;
		int horizontalScore=cellToLeft.getScore()-gapPenalty;
		int diagonalScore=cellAboveLeft.getScore();
		diagonalScore=diagonalScore+matchScore;
		if(verticalScore>=horizontalScore) {
			if(verticalScore>=diagonalScore) {
				currentCell.setScore(verticalScore);
				currentCell.setPrevCell(cellAbove);
			}
			else {
				currentCell.setScore(diagonalScore);
				currentCell.setPrevCell(cellAboveLeft);
			}
		}
		else {
			if(horizontalScore>=diagonalScore) {
				currentCell.setScore(horizontalScore);
				currentCell.setPrevCell(cellToLeft);
			}
			else {
				currentCell.setScore(diagonalScore);
				currentCell.setPrevCell(cellAboveLeft);
			}
		}
	}
	/*
	 * protected static Cell[][] fillInTable(String str1, String str2,int choice): 
	 * This method fills in the whole alignment matrix depending on the input choice on DNA/protein sequence alignment. 
	 */	
	protected static Cell[][] fillInTable(String str1, String str2,int choice) {
	    	Cell[][] scoreTable=new Cell[str2.length()+1][str1.length()+1];
	    	initialize(scoreTable);
	    	for(int i=1;i<=str2.length();i++) {
	    		for(int j=1;j<=str1.length();j++) {
	    			if(choice==0) {
	    			fillInCellDNASeq(str1,str2, scoreTable[i][j],scoreTable[i-1][j], scoreTable[i][j-1], scoreTable[i-1][j-1]);
	    			}
	    			else {
	    				fillInCellProSeq(str1,str2, scoreTable[i][j],scoreTable[i-1][j], scoreTable[i][j-1], scoreTable[i-1][j-1]);
	    			}
	    		}
	    	}
	    	return scoreTable;
	}
	/*
	 * protected static String[] trackBack(String str1, String str2,Cell[][] scoreTable): 
	 * This method constructs the alignment from tracing back path. 
	 * Tracing back begins at the element in the bottom right hand corner of the matrix and 
	 * then it follows pointers that gave maximum score for each element. 
	 * Tracing back continues until the top left hand corner of the matrix is reached.
	 */
	protected static String[] trackBack(String str1, String str2,Cell[][] scoreTable) {
		StringBuffer Buf1=new StringBuffer();
		StringBuffer Buf2=new StringBuffer();
		Cell startingCell=scoreTable[scoreTable.length-1][scoreTable[0].length-1];
		//Score of alignment will be value in element bottom right hand corner
		System.out.println("Score of alignment = "+startingCell.getScore());
		Cell currentCell=startingCell;
		while(currentCell.getPrevCell() != null) { //   currentCell.getCol()!=0 && currentCell.getRow()!=0 CANNOT BE USED HERE!
			if(currentCell.getRow()-currentCell.getPrevCell().getRow()==1) {
				Buf2.insert(0, str2.charAt(currentCell.getRow()-1));
			}
			else {
				Buf2.insert(0, '-');
			}
			if(currentCell.getCol()-currentCell.getPrevCell().getCol()==1) {
				Buf1.insert(0, str1.charAt(currentCell.getCol()-1));
			}
			else {
				Buf1.insert(0, '-');
			}
//			if(currentCell.getRow()-currentCell.getPrevCell().getRow()==1 && currentCell.getCol()-currentCell.getPrevCell().getCol()!=1) {//up arrow
//				Buf1.insert(0, '-');
//				Buf2.insert(0, str2.charAt(currentCell.getRow()-1));
//			}
//			if(currentCell.getCol()-currentCell.getPrevCell().getCol()==1 && currentCell.getRow()-currentCell.getPrevCell().getRow()!=1) {//left arrow
//				Buf2.insert(0, '-');
//				Buf1.insert(0, str1.charAt(currentCell.getCol()-1));
//			}
//			if(currentCell.getCol()-currentCell.getPrevCell().getCol()==1 && currentCell.getRow()-currentCell.getPrevCell().getRow()==1) {
//				Buf1.insert(0, str1.charAt(currentCell.getCol()-1));
//				Buf2.insert(0, str2.charAt(currentCell.getRow()-1));
//			}
			
			currentCell=currentCell.getPrevCell();
		}
		String[] alignments=new String[] {Buf1.toString(), Buf2.toString()};
		return alignments;
	}
	/*
	 * protected static void repeatedProcess(String strt1, String strt2,int choice): This method represents the process of aligning two sequences
	 *  based on the input choice parameter. It first fills in the alignment matrix and then traces back to get the final pairwise alignment.
	 *  When choice=0, it aligns the DNA sequences; when choice=1, it aligns the protein sequences. 
	 */
	protected static void repeatedProcess(String strt1, String strt2,int choice) {
		System.out.println("sequence 1:"+strt1);
		System.out.println("sequence 2:"+strt2);
		Cell[][] scoreTablet1=fillInTable(strt1, strt2,choice);	
		System.out.println("Below is the global seguence alignment table:");
		for(int j=0;j<=strt2.length();j++) {
			for(int i=0;i<=strt1.length();i++) {
				System.out.printf("%5d",scoreTablet1[j][i].getScore());
			}
			System.out.println();
		}
		String[] alignmentst1=trackBack(strt1, strt2, scoreTablet1);
		System.out.println("Below is the global seguence alignment result:");
		System.out.println("Sequence 1");
		System.out.println(alignmentst1[0]);
		System.out.println("Sequence 2");
		System.out.println(alignmentst1[1]);
	}
	/*
	 * protected static void TestingMode(): This method runs three testing cases respectively for DNA sequence alignment 
	 * and protein sequence alignment. 
	 */
	protected static void TestingMode() {
		System.out.println();
		System.out.println();
		System.out.println("************************************************************************");
		System.out.println("Testing mode begins!");
		System.out.println("DNA sequence alignment testing:");
		System.out.println("-------------Test 1 Begins--------------");
		String strt1="ACCGTA";
		String strt2="ACGTT";
		repeatedProcess(strt1,strt2,0);
		System.out.println("--------------Test 1 Ends----------------");
		System.out.println("-------------Test 2 Begins--------------");
		String strt3="CGGATTGCATTACG";
		String strt4="CGATTCATG";
		repeatedProcess(strt3,strt4,0);
		System.out.println("--------------Test 2 Ends----------------");
		System.out.println("-------------Test 3 Begins--------------");
		String strt5="TAGGCTGAGCGACCGCGTA";
		String strt6="TGGCCTAGGCCGA";
		repeatedProcess(strt5,strt6,0);
		System.out.println("--------------Test 3 Ends----------------");
		
		System.out.println("Protein sequence alignment testing:");
		System.out.println("-------------Test 1 Begins--------------");
		String str1="MEKVNEERDAVF";
		String str2="EDHIGDRRRSV";
		repeatedProcess(str1,str2,1);
		System.out.println("--------------Test 1 Ends----------------");
		System.out.println("-------------Test 2 Begins--------------");
		String str3="RSLLEEA";
		String str4="FADEMEKT";
		repeatedProcess(str3,str4,1);
		System.out.println("--------------Test 2 Ends----------------");
		System.out.println("-------------Test 3 Begins--------------");
		String str5="SYDVEVADTPQPHIPIRFRHPPIA";
		String str6="MKRGQAVDFCHWVSHLIATEIDEK";
		repeatedProcess(str5,str6,1);
		System.out.println("--------------Test 3 Ends----------------");
		
		System.out.println("Testing mode ends and follow the guidance you can try whatever you want!");
		System.out.println("************************************************************************");
		System.out.println();
		System.out.println();
	}
	
	/*
	 * public static void main(String[] args) throws IOException: Firstly, the main function runs the six testing cases. 
	 * Afterwards, the user can test whatever they like based on their preference on whether aligning DNA sequences or 
	 * protein sequences by changing the value for the parameter choice. 
	 */
	public static void main(String[] args) throws IOException {
		int continueOrNot=0;
		System.out.println("***************************************************");
		System.out.println("Welcome to Jin Tao's Sequence Analysis Application!");
		System.out.println("*              Student ID 11474660                *");
		System.out.println("***************************************************");
		//testing mode
		TestingMode();
		//ask for choice
		while(continueOrNot==0) {
		Scanner scanner = new Scanner(System.in);
	    System.out.println("What kind of sequence would you like to align? Enter 0 for DNA sequence; 1 for protein sequence: ");
	    int choice = scanner.nextInt();
	    System.out.printf("You choose " + choice);
	    if(choice==0) {
	    	System.out.println(" to align DNA sequences!");
	    }
	    else {
	    	System.out.println(" to align protein sequences!");
	    }
		BufferedReader br=new BufferedReader(new InputStreamReader(System.in));
		System.out.println("Enter the sequence 1:");
		String str1=br.readLine();
		System.out.println("Enter the sequence 2:");
		String str2=br.readLine();
		repeatedProcess(str1, str2,choice);
		
		System.out.println("Do you want to continue aligning sequences? Enter 0 to continue, or enter 1 to exit.");
		Scanner sc = new Scanner(System.in);
		continueOrNot = sc.nextInt();
		}
		System.out.println("Goodbye!");
	}
}
	


