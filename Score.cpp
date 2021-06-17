#include "Score.h"

/*

 */
double Score::shScore(vector<int> &mismatches, vector<string> &hsuKeys, map<string, vector<double>> &hsuMatrix)
{
	double tot_sh = 1.0;
	for (int i = 0; i < mismatches.size(); i++)
	{
		tot_sh *= hsuMatrix[hsuKeys[i]][mismatches[i]];
	}
	return tot_sh;
}

/*

 */
double Score::ssScore(vector<int> &mismatches)
{
	double tot_ss = 1.0;
	for (int i = 0; i < mismatches.size(); i++)
	{
		if (mismatches[i] < 6)
		{
			tot_ss -= 0.1;
		}
		else if (mismatches[i] < 12)
		{
			tot_ss -= 0.05;
		}
		else
		{
			tot_ss -= 0.0125;
		}
	}
	return tot_ss;
}

/*

 */
double Score::stScore(vector<int> &mismatches)
{
	double tot_st = 3.5477;
	for (int i = 0; i < mismatches.size(); i++)
	{
		tot_st -= 1.0 / (mismatches[i] + 1);
	}
	return tot_st / 3.5477;
}