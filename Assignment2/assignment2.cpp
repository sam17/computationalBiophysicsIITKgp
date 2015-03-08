#include <stdio.h>
#include <string>
#include <algorithm>
#include <iostream>
#include <vector>

using namespace std;

#define MAX(a,b) (a)>(b)?(a):(b)
//BLOSUM62 Matrix
static int blosum[20][20] = {
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
	{ 0, -3, -3, -3, -1, -2, -2, -3, -3,  3,  1, -2,  1, -1, -2, -2,  0, -3, -1,  4}};

//Function to get ID from Protein sequence character for BLOSUM62 matrix 
int getID(char ch)
{
    switch(ch)
    {
    case 'A': return 0;
    case 'R': return 1;
    case 'N': return 2; 
    case 'D': return 3;
    case 'C': return 4;
    case 'Q': return 5;
    case 'E': return 6;
    case 'G': return 7;
    case 'H': return 8;
    case 'I': return 9;
    case 'L': return 10;
    case 'K': return 11;
    case 'M': return 12;
    case 'F': return 13;
    case 'P': return 14;
    case 'S': return 15;
    case 'T': return 16;
    case 'W': return 17;
    case 'Y': return 18;
    case 'V': return 19;
    default: return -1;
    }
}

//Function to get the match.mismatch score from BLOSUM62 Matrix
int getBLOSUMscore(char a,char b)
{
    int id_a = getID(a);
    int id_b = getID(b);
    
    return blosum[id_a][id_b];
}

int main(int argc,char *argv[])
{
    string sequenceA,sequenceB;
    double gap_weight;

    //Input Sequence and Gap Weight
    cin>>sequenceA;
    cin>>sequenceB;
    cin>>gap_weight;

    //DP Matrix
    vector<vector<int> > matrix(sequenceB.size()+1, std::vector<int>(sequenceA.size()+1));

    //Initializations
    for(int i=0;i<=sequenceB.size();i++)
	matrix[i][0] = gap_weight*i;
    for(int i=0;i<=sequenceA.size();i++)
	matrix[0][i] = gap_weight*i;

    //Matrix Propagation
    for(int i=1;i<=sequenceB.size();i++)
    {
	for(int j=1;j<=sequenceA.size();j++)
	{
	    int a = matrix[i-1][j-1]+getBLOSUMscore(sequenceB[i-1],sequenceA[j-1]);
	    int b = matrix[i][j-1]+gap_weight;
	    int c = matrix[i-1][j]+gap_weight;
	    matrix[i][j] = MAX(a,MAX(b,c));
	}
    }
    
    /*
    for(int i=0;i<=sequenceB.size();i++)
    {
	for(int j=0;j<=sequenceA.size();j++)
	{
	    cout<<matrix[i][j]<<" ";
	}
	cout<<endl;
    }
    */

    //Backtracing to get aligned sequences
    string AlignmentA = "";
    string AlignmentB = "";
    int j = sequenceA.size();
    int i = sequenceB.size();
    
    cout<<sequenceA<<endl<<sequenceB<<endl<<"Alignment Score: "<<matrix[i][j]<<endl;
    
    while(i>0 || j>0)
    {
	if (i > 0 && j > 0 && (matrix[i][j] == matrix[i-1][j-1] + getBLOSUMscore(sequenceA[j-1], sequenceB[i-1])))
	{
	    AlignmentA = sequenceA[j-1] + AlignmentA;
	    AlignmentB = sequenceB[i-1] + AlignmentB;
	    i = i - 1;
	    j = j - 1;
	}
	else if (i > 0 && matrix[i][j] == (matrix[i-1][j] + gap_weight))
	{
	    AlignmentB = sequenceB[i-1] + AlignmentB;
	    AlignmentA = "-" + AlignmentA;
	    i = i - 1;
	}
	else if (j > 0 and matrix[i][j] == (matrix[i][j-1] + gap_weight))
	{
	    AlignmentB = "-" + AlignmentB;
	    AlignmentA = sequenceA[j-1] + AlignmentA;
	    j = j - 1;
	}
	else cout<<"ERR\n";
    }
    
    cout<<AlignmentA<<endl<<AlignmentB;
    
    return 0;
}
