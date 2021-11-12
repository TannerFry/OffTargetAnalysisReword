#include "Score.h"

#include <iostream>
using namespace std;

/*
	function for calculating the sh score

	@param mismatches	=> locations of the mismatches found between 2 sequences
	@param hsuKeys		=> keys for indexing into the HSU matrix
	@param hsuMatrix	=> 2D HSU matrix
 	@param gRNA_length	=> length of gRNA sequence

	@return tot_sh	=> final sh score

*/
double Score::shScore(vector<int> &mismatches, vector<string> &hsuKeys, map<string, vector<double>> &hsuMatrix, int &gRNA_length)
{
	double tot_sh = 1.0;
	for (int i = 0; i < mismatches.size(); i++)
	{
		tot_sh *= hsuMatrix[hsuKeys[i]][gRNA_length - mismatches[i]];
	}
	return tot_sh;
}

/*
	function for calculating the ss score

	@param mismatches	=> locations of the mismatches found between 2 sequences
	@param gRNA_length	=> length of gRNA sequence

	@return tot_ss	=> final ss score

*/
double Score::ssScore(vector<int> &mismatches, int &gRNA_length)
{
	double tot_ss = 1.0;
	
	if (gRNA_length == 24)
	{
		for (int i = 0; i < mismatches.size(); i++)
		{
			if (mismatches[i] <= 8)
			{
				tot_ss -= 0.1;
			}
			else if (mismatches[i] <= 20)
			{
				tot_ss -= 0.0125;
			}
			else
			{
				tot_ss -= 0;
			}
		}
	}
	else
	{
		for (int i = 0; i < mismatches.size(); i++)
		{
			if (mismatches[i] <= 6)
			{
				tot_ss -= 0.1;
			}
			else if (mismatches[i] <= 12)
			{
				tot_ss -= 0.05;
			}
			else
			{
				tot_ss -= 0.0125;
			}
		}
	}

	return tot_ss;
}

/*
	function for calculating the st score

	@param mismatches	=> locations of the mismatches found between 2 sequences

	@return tot_st	=> final st score

*/
double Score::stScore(vector<int> &mismatches)
{
	double tot_st = 3.547;
	for (int i = 0; i < mismatches.size(); i++)
	{
		tot_st -= (1.0 / (mismatches[i]));
	}
	return tot_st / 3.5477;
}