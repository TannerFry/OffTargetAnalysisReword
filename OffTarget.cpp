#include "OffTarget.h"

/* 
	function for parsing input arguments 

	@param argc	=> number of input arguments
	@param argv	=> array of input arguments
*/
void OffTarget::parseInputArguments(int argc, char *argv[])
{
	/* parse and store input args */
	queryFilePath = string(argv[1]);
	endo = string(argv[2]);
	csprFilePath = string(argv[3]);
	sqlFilePath = string(argv[4]);
	outputFilePath = string(argv[5]);
	casperInfoFilePath = string(argv[6]);
	maxMismatches = stoi(string(argv[7]));
	threshold = stod(string(argv[8]));
	if (string(argv[9]).find('T') != string::npos) {
		detailedOutput = true;
	}
	if (string(argv[10]).find('T') != string::npos) {
		avgOutput = true;
	}
	hsuMatrixName = string(argv[11]);
}

/*
	function for calling file operations object to parse data needed for algorithm:
*/
void OffTarget::getAlgorithmData()
{
	FileOp.parseCasperInfo(casperInfoFilePath, endo, endoData, hsuMatrixName, hsuMatrix);
	FileOp.parseCsprFile(csprFilePath, uniqueSeqs, uniqueScores, uniqueLocations, uniqueChroms);
	FileOp.parseSqlFile(sqlFilePath, repeatSeqs, repeatScores, repeatLocations, repeatChroms);
	FileOp.parseQueryFile(queryFilePath, querySeqs, queryScores);
}

/* 
	function for running the OffTarget algorithm
 */
void OffTarget::run()
{
	/* vars */
	int seqLength = endoData[4];
	string currentQuerySeq = "", details = "\n";
	vector<double> targets;
	int threadCount = queryScores.size();
	vector<thread> runningThreads(threadCount);
	int currentQueryScore = 0;
	
	/* init multi-threading variables */
	vector<vector<vector<double> > > targetScores(queryScores.size(), vector<vector<double> >(2));
	vector<vector<vector<unsigned long > > > targetIndexes( queryScores.size(), vector<vector<unsigned long> >(2));
	
	/* open output file object so findSimilars can write out results*/
	FileOp.openOutputFile(outputFilePath, avgOutput);
	/* score each query sequence and write out the results in a thread */
	for (unsigned long i = 0; i < queryScores.size(); i++)
	{
		/* get current query variables for thread */
		currentQuerySeq = querySeqs.substr(i * seqLength, seqLength);
		currentQueryScore = queryScores[i];
		
		/* create thread */
		thread t([this, currentQuerySeq, currentQueryScore, &targetScores, &targetIndexes, i]()
		{ 
			findSimilars(currentQuerySeq, currentQueryScore, targetScores[i], targetIndexes[i]); 
		});

		/* save thread */
		runningThreads[i] = move(t);
	}

	/* join threads and write out data to file */
	for (unsigned long i = 0; i < runningThreads.size(); i++)
	{
		runningThreads[i].join();
		currentQuerySeq = querySeqs.substr(i * seqLength, seqLength);

		/* write findings */
		FileOp.writeResults(avgOutput, currentQuerySeq, targetScores[i][0], targetIndexes[i][0], targetScores[i][1], targetIndexes[i][1], repeatLocations, repeatChroms, repeatSeqs, uniqueLocations, uniqueChroms, uniqueSeqs, seqLength);
		
		/* clear out data */
		targetScores[i].clear();
		targetScores[i].shrink_to_fit();
		targetIndexes[i].clear();
		targetIndexes[i].shrink_to_fit();
	}

	/* close output file */
	FileOp.closeOutputFile();
}

/*
	OffTarget function for finding similar sequences in the reference organism, scoring the findings, and writing out the results

	@param currentQuerySeq		=> current query sequence being analyzed
	@param currentQueryScore	=> on-target score of query sequence
 */
void OffTarget::findSimilars(string currentQuerySeq, int currentQueryScore, vector<vector<double> > &targetScores, vector<vector<unsigned long> > &targetIndexes)
{
	//vars
	int seqLength = endoData[4];

	/* run query sequence against unique sequences from CSPR file */
	findSimilarsUnique(currentQuerySeq, currentQueryScore, seqLength, targetScores[0], targetIndexes[0]);

	/* run query sequence against repeat sequences from DB file */
	findSimilarsRepeat(currentQuerySeq, currentQueryScore, seqLength, targetScores[1], targetIndexes[1]);
}

/*
	function for running off target analysis of query sequence against the unique organism data from CSPR file

	@param currentQuerySeq		=> current query sequence string
	@param currentQueryScore	=> on-score of current query string
	@param seqLength			=> length of sequences for current endo
	@param targetScores			=> vector that gets filled with individual scores from targets found in this function
	@param uniqueIndexes		=> vector that gets filled with the index values of target sequences found

*/
void OffTarget::findSimilarsUnique(string &currentQuerySeq, int &currentQueryScore, int &seqLength, vector<double> &targetScores, vector<unsigned long> &targetIndexes)
{	
	/* vars */
	double rRatio = 0.0;
	double value = 0.0;
	string refSeq;

	/* loop through each organims unique seq from CSPR file and compare against the given query sequence */
	for (unsigned long i = 0; i < uniqueSeqs.length() / seqLength; i++)
	{
		vector<int> mismatches;
		vector<string> hsuKeys;
		rRatio = uniqueScores[i] / currentQueryScore;
		refSeq = uniqueSeqs.substr(i * seqLength, seqLength);

		//character by character comparison of ref and query sequences
		if (getMismatches(refSeq, currentQuerySeq, mismatches, hsuKeys, seqLength))
		{
			//update runnning score if mismatch count wasn't too large
			targetIndexes.push_back(i);
			value = (((sqrt(score.shScore(mismatches, hsuKeys, hsuMatrix)) + score.stScore(mismatches)) * pow(score.ssScore(mismatches), 6) * pow(rRatio, 2)) / 4);
			targetScores.push_back(value);
		}
	}
}

/* 
	function for running off target analysis of query sequence against the repeat organism data from DB file
*/
void OffTarget::findSimilarsRepeat(string &currentQuerySeq, int &currentQueryScore, int &seqLength, vector<double> &targetScores, vector<unsigned long> &targetIndexes)
{
	/* vars */
	double rRatio = 0.0;
	double value = 0.0;
	string refSeq;

	/* loop through each organims unique seq from CSPR file and compare against the given query sequence */
	for (unsigned long i = 0; i < repeatSeqs.length() / seqLength; i++)
	{
		vector<int> mismatches;
		vector<string> hsuKeys;
		rRatio = repeatScores[i] / currentQueryScore;
		refSeq = repeatSeqs.substr(i * seqLength, seqLength);

		//character by character comparison of ref and query sequences
		if (getMismatches(refSeq, currentQuerySeq, mismatches, hsuKeys, seqLength))
		{
			//update runnning score if mismatch count wasn't too large
			targetIndexes.push_back(i);
			value = (((sqrt(score.shScore(mismatches, hsuKeys, hsuMatrix)) + score.stScore(mismatches)) * pow(score.ssScore(mismatches), 6) * pow(rRatio, 2)) / 4);
			targetScores.push_back(value);
		}
	}
}

/*
	function for character wise comparison of two sequences

	@param refSeq				=> sequence pulled from organism CSPR/DB file
	@param currentQuerySeq		=> sequence pulled from query file
	@param mismatchLocations	=> vector that will get filled in with the locations the mismatched characters appear at between the two sequences
	@param hsuKeys				=> vector that will get filled in with the HSU keys to use with each mismatch location
	@param seqLength			=> int representing the length of sequences for the current endo

	@return true	=> number of mismatches found not larger than max mismatches allowed
	@return false	=> otherwise
 */
bool OffTarget::getMismatches(string &refSeq, string &currentQuerySeq, vector<int> &mismatchLocations, vector<string> &hsuKeys, int &seqLength)
{
	/* vars */
	string hsuKey;
	
	//character by character comparison of ref and query sequences
	for (int j = seqLength - 1; j >= 0; j--)
	{
		try
		{
			if (refSeq.at(j) != currentQuerySeq.at(j))
			{
				//store mismatch location
				mismatchLocations.push_back(j);

				//store key for HSU matrix
				hsuKey = string() + currentQuerySeq.at(j) + reverseComp(refSeq.at(j));
				hsuKeys.push_back(hsuKey);
			}
		}
		catch (int e)
		{
			cout << "mismatch error" << endl;
		}
		
		/* if there are too many mismatches, break */
		if (mismatchLocations.size() > maxMismatches)
		{
			return false;
		}
	}
	return true;
}

/*
	function to reverse complement a character
	
	@param c => character to perfrom reverse complement on 
	
	@return reverse complement of input character
 */
char OffTarget::reverseComp(char &c)
{
	char n = 'N';
	switch (c)
	{
	case 'A':
		n = 'T';
		break;
	case 'T':
		n = 'A';
		break;
	case 'G':
		n = 'C';
		break;
	case 'C':
		n = 'G';
		break;
	}
	return n;
}