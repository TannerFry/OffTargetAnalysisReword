#pragma once
#include <vector>
#include <cstdint>
#include <string>
#include <map>

using namespace std;

class Score
{
	public:
		/* function for calculating the sh score */
		double shScore(vector<int> &mismatches, vector<string> &hsuKeys, map<string, vector<double >> &hsuMatrix, int &gRNA_length);

		/* function for calculating the ss score */
		double ssScore(vector<int> &mismatches, int &gRNA_length);

		/* function for calculating the st score */
		double stScore(vector<int> &mismatches);
};
