#pragma once
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <sstream>
#include <cstdint>
#include <numeric>
#include <iomanip>
#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include "sqlite3.h"

using namespace std;

class FileOperations
{
	public:
		/* functin for parsing CASPERinfo file to retrieve endo data and HSU matrix */
		void parseCasperInfo(string &casperInfoFile, string &endo, vector<int> &endoData, string &hsuMatrixName, map<string, vector<double>> &hsuMatrix);
		
		/* function for parsing organism CSPR file to retrieve the unique reference targets */
		void parseCsprFile(string &csprFilePath, string &uniqueSeqs, vector<uint8_t> &uniqueScores, vector<long long> &uniqueLocations, vector<int> &uniqueChroms);

		/* function for parsing organism SQL file to retrieve the repeat reference targets */
		void parseSqlFile(string &dbFilePath, string &repeatSeqs, vector<uint8_t> &repeatScores, vector<long long> &repeatLocations, vector<int> &repeatChroms);
		
		/* function for parsing the input query sequences to score */
		void parseQueryFile(string &queryFilePath, string &querySeqs, vector<uint8_t> &queryScores);
		
		/* function to open output file */
		void openOutputFile(string &outputFilePath, bool &avgOutput);

		/* function to close output file */
		void closeOutputFile();

		/* function to write out the scoring results */
		void writeResults(bool &avgOutput, string &querySeq, vector<double> &uniqueScores, vector<unsigned long> &uniqueIndexes, vector<double> &repeatScores, vector<unsigned long> &repeatIndexes, vector<long long> &repeatLocations, vector<int> &repeatChroms, string &repeatSeqs, vector<long long> &uniqueLocations, vector<int> &uniqueChroms, string &uniqueSeqs, int &seqLength);
	
	private:
		ofstream outputFile;


		/* hsuKeys => static keys for the HSU matrix */
		vector<string> hsuKeys = {"GT", "AC", "GG", "TG", "TT", "CA", "CT", "GA", "AA", "AG", "TC", "CC"};
		
		/* function split a given string based on a delimter */
		vector<string> split(string s, char &delimiter);

		/* function to load the default HSU MATRIX-spCas9-2013 */
		void loadDefaultHsuMatrix(map<string, vector<double>> &hsuMatrix);
};
