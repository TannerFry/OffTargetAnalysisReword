#pragma once
#include "FileOperations.h"
#include "Score.h"
#include <thread>

using namespace std;

/* OffTarget class represents the primary object of the algorithms implementation */
class OffTarget
{
	public:
		/* function for parsing input arguments */
		void parseInputArguments(int argc, char *argv[]);
		
		/* parse data needed for algorithm */
		void getAlgorithmData();
		
		/* function for running the OffTarget algorithm */
		void run();
	
	private:
		/*	
			Input argument variable definitions:
			queryFile		=> File path to input file holding sequence data to run against the reference organism data and get off target score for
			endo			=> Defines endonuclease used 
			csprFile		=> File path to GZIP cspr file of organism
			sqlFile			=> File path for SQL repeats file of organism
			outputFile		=> File path for output file
			casperInfoFile	=> File path to CASPERinfo
			maxMismatches	=> Defines the number of max mismatched letters between two sequences
			threshold		=> Defines the threshold value for scores: score below this threshold will not be reported
			avgOutput		=> Defines if average output format is used: only the offtarget score is shown for each query sequences in the output file
			detailedOutput	=> Defines if detailed output format is used: provides additional information on targets found from the algorithm
			hsuMatrixName	=> HSU matrix name to parse from CASPERinfo
			three_prime		=> boolean - True = 3 prime, False = 5 prime
		*/
		bool avgOutput = false, detailedOutput = false;
		int maxMismatches = 0;
		double threshold = 0;
		string endo, queryFilePath, csprFilePath, sqlFilePath, outputFilePath, casperInfoFilePath, hsuMatrixName;
		bool three_prime = true;

		/* 
			CASPERinfo variable definitions
			endoData	=> vector containing {pam length, 3' length, seed length, 5' length, sequence length}
			hsuMatrix	=> See CASPERinfo for HSU matrix structure
		*/
		vector<int> endoData;
		map<string, vector<double>> hsuMatrix;
		vector<string> hsuKeys = { "GT", "AC", "GG", "TG", "TT", "CA", "CT", "GA", "AA", "AG", "TC", "CC" };

		/*
			CSPR file variable definitions
			uniqueSeqs		=> concatenated string of all unqiue sequences in CSPR file
			uniqueScores	=> vector of ints holding the scores of each unqiue sequence
			uniqueLocations	=> vector holding the locations of sequences from CSPR file
			uniqueChroms	=> vector holding the chromosome of each sequence from CSPR file
		*/
		string uniqueSeqs;
		vector<uint8_t> uniqueScores;
		vector<long long> uniqueLocations;
		vector<int> uniqueChroms;

		/*
			DB file variable definitions

		*/
		string repeatSeqs;
		vector<uint8_t> repeatScores;
		vector<long long> repeatLocations;
		vector<int> repeatChroms;

		/*
			Query file variable definitions
			querySeqs	=> concatenated string of all sequences in the query file
			queryScores	=> vector of ints holding the scores of each query sequence
		*/
		string querySeqs;
		vector<uint8_t> queryScores;

		/* FileOperations object - used for all file parsing/writing operations */
		FileOperations FileOp;

		/* score object to run scoring algorithms */
		Score score;

		/* vector to hold OffTarget scores for query sequences */
		vector<double> queryOffTargetScores;

		/* 	OffTarget analysis function for finding similar sequences in the reference organism, scoring the findings, and writing out the results 
			findSimilars is a wrapper for calling findSimilarsUnique and findSimiarsRepeat for each query sequence
		*/
		void findSimilars(string currentQuerySeq, int currentQueryScore, vector<vector<double> > &targetScores, vector<vector<unsigned long> > &targetIndexes);

		/* function for running off target analysis of query sequence against the unique organism data from CSPR file */
		void findSimilarsUnique(string &currentQuerySeq, int &currentQueryScore, int &seqLength, vector<double> &runningScores, vector<unsigned long> &targetIndexes);

		/* function for running off target analysis of query sequence against the repeat organism data from DB file */
		void findSimilarsRepeat(string &currentQuerySeq, int &currentQueryScore, int &seqLength, vector<double> &runningScores, vector<unsigned long> &targetIndexes);

		/* function for character wise comparison of two sequences */
		bool getMismatches(string &refSeq, string &currentQuerySeq, vector<int> &mismatchLocations, vector<string> &hsuKeys, int &seqLength);

		/* function to reverse complement a character */
		char reverseComp(char &c);
};
