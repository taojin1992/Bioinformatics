/*
***************************************************************************
*Jin Tao
*Student ID 11474660
Class Description:
Cell Class:
Each object in the Cell class represents one element in the alignment matrix. 
According to the Global Sequence Alignment Dynamic Programming Algorithm, 
each element has a score and a pointer for final tracking back. 
*****************************************************************************
*/
public class Cell {
	private int row;//Instance variable row and col are used to specify the position of each Cell object in the alignment matrix. 
	private int col;
	private Cell prevCell;//Instance variable prevCell represents the tracking back pointer of each Cell object. 
	private int score;//Instance variable score represents the score of each Cell object. 
	/*
	 * public Cell(int i, int j) creates a Cell object according to the position in the alignment matrix, 
	 * with 0 as the initial score and null as the initial tracking back pointer.
	 */
	public Cell(int i, int j)
	{
	    prevCell=null;
		row = i;
	    col = j;
	    score=0;
	}
	//public void setScore(int Score): This setter method sets the score for each Cell object. 
	public void setScore(int Score) {
		score=Score;		
	}
	//public void setPrevCell(Cell Pointer): This setter method sets the arrow inside each Cell object for back tracking. 
	public void setPrevCell(Cell Pointer) {
		prevCell=Pointer;
	}
	//public int getScore(): This getter method gets the score of each Cell object.
	public int getScore() {
		return score;
	}
	//public int getRow(): This getter method gets the row number of each Cell object in the alignment matrix.
	public int getRow() {
		return row;
	}
	//public int getCol(): This getter method gets the column number of each Cell object in the alignment matrix.
	public int getCol() {
		return col;
	}
	//public Cell getPrevCell():This getter method gets the arrow inside each Cell object in the alignment matrix. 
	public Cell getPrevCell() {
		return prevCell;
	}	
}
