#pragma once
#include <vector>
#include <cstdint>
#include <string>
#include <map>

using namespace std;

class Score
{
	public:
		/**/
		double shScore(vector<int> &mismatches, vector<string> &hsuKeys, map<string, vector<double >> &hsuMatrix, int &gRNA_length);

		/**/
		double ssScore(vector<int> &mismatches, int &gRNA_length);

		/**/
		double stScore(vector<int> &mismatches);
};
