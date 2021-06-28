#include "FileOperations.h"

/*
	functin for parsing CASPERinfo file to retrieve endo data and HSU matrix

	@param casperInfoFilePath	=> file path to CASPERinfo
	@param endo					=> name of endo being used
	@param endoData				=> vector to hold endo data for OffTarget class
	@param hsuMatrixName		=> name of HSU matrix to be loaded
	@param hsuMatrix			=> map to hold the HSU matrix data for OffTarget class

 */
void FileOperations::parseCasperInfo(string &casperInfoFilePath, string &endo, vector<int> &endoData, string &hsuMatrixName, map<string, vector<double>> &hsuMatrix)
{
	/* vars */
	ifstream casperInfoFile;
	string line;
	char endoDelimeter = ';';
	double hsuVal = 0;
	bool foundEndo, foundHsuMatrix = false;

	//open and confirm CASPERinfo opened
	casperInfoFile.open(casperInfoFilePath);
	if (casperInfoFile.is_open())
	{
		//parse CASPERinfo for the endo data
		while (getline(casperInfoFile, line))
		{
			if (line.find(endo + ";") != string::npos)
			{
				vector<string> splitLine = split(line, endoDelimeter);
				endoData.push_back(splitLine[1].length()); //pam length
				endoData.push_back(stoi(splitLine[2]));	//3' length
				endoData.push_back(stoi(splitLine[3])); //seed length
				endoData.push_back(stoi(splitLine[4])); //5' length
				endoData.push_back(stoi(splitLine[2]) + stoi(splitLine[3]) + stoi(splitLine[4])); //sequence length
				foundEndo = true;
				break;
			}
		}

		//parse CASPERinfo for the HSU matrix data
		casperInfoFile.seekg(SEEK_SET);
		while (getline(casperInfoFile, line))
		{
			if (line.find(hsuMatrixName) != string::npos)
			{
				for (int i = 0; i < 12; i++)
				{
					getline(casperInfoFile, line);
					stringstream ss(line);
					while (ss >> hsuVal)
					{
						hsuMatrix[hsuKeys[i]].push_back(hsuVal);
					}
				}
				foundHsuMatrix = true;
				break;
			}
		}
		casperInfoFile.close();
	}
	/* exit if CASPERinfo could not be found */
	else
	{
		cerr << "CASPERinfo file could not be opened." << endl;
		exit(-1);
	}

	/* exit if endo data couldn't be found in CASPERinfo */
	if (foundEndo == false)
	{
		cerr << "Endo information not found in CASPERinfo." << endl;
		exit(-1);
	}

	/* if HSU matrix specified couldn't be found, load the default setup */
	if (foundHsuMatrix == false)
	{
		loadDefaultHsuMatrix(hsuMatrix);
	}
}

/*
	function for parsing organism CSPR file to retrieve the unique reference targets

	@parma csprFilePath	=> file path to CSPR file
	@param uniqueSeqs	=> concatenated string of all unqiue sequences in CSPR file
	@param uniqueScores	=> vector of ints holding the scores of each unqiue sequence
 */
void FileOperations::parseCsprFile(string &csprFilePath, string &uniqueSeqs, vector<uint8_t> &uniqueScores, vector <long long > &uniqueLocations, vector<int> &uniqueChroms)
{
	/* vars */
	string line;
	char csprFileDelimeter = ',';
	int chromCount = 0;

	/* open and verify CSPR file */
	ifstream file(csprFilePath, ios_base::in | ios_base::binary);
	if (file.is_open())
	{
		boost::iostreams::filtering_streambuf<boost::iostreams::input> inbuf;
		inbuf.push(boost::iostreams::gzip_decompressor());
		inbuf.push(file);
		istream instream(&inbuf);

		/* ignore CSPR file meta data */
		getline(instream, line);
		getline(instream, line);
		getline(instream, line);

		/* loop through each line and extract the sequence data */
		while (getline(instream, line))
		{
			if (line.find(">") == string::npos)
			{
				vector<string> splitLine = split(line, csprFileDelimeter);
				uniqueSeqs += splitLine[1];
				uniqueScores.push_back(stoi(splitLine[3]));
				uniqueLocations.push_back(stoll(splitLine[0]));
				uniqueChroms.push_back(chromCount);
			}
			else
			{
				chromCount++;
			}
		}
	}
	else
	{
		cerr << "CSPR file could not be opened." << endl;
		exit(-1);
	}
}


/*
	function for parsing organism SQL file to retrieve the repeat reference targets
 */
void FileOperations::parseSqlFile(string &dbFilePath, string &repeatSeqs, vector<uint8_t> &repeatScores, vector<long long> &repeatLocations, vector<int> &repeatChroms)
{
	// the sql we will turn in to a prepared statement
	string sql = "SELECT seed, chromosome, location, three, five, score FROM repeats;";
	sqlite3 *db; // pointer to our db
	sqlite3_stmt *pstmt; // prepared statements corresponding to sql
	string seed, chromosome, location, three, five, score;
	char delimeter = ',';
	char *zErrMsg = 0;
	int rc; // return code from sqlite library

	// open the database
	rc = sqlite3_open_v2(dbFilePath.c_str(), &db, SQLITE_OPEN_READONLY, NULL);
	if (rc)
	{
		fprintf(stderr, "Can't open database: %s\n", sqlite3_errmsg(db));
		sqlite3_close_v2(db);
		exit(-1);
	}

	// create a prepared statement
	rc = sqlite3_prepare_v3(db, sql.c_str(), -1, 0, &pstmt, NULL);
	if (rc)
	{
		fprintf(stderr, "Couldn't prepare sql statement: %s\n", sqlite3_errmsg(db));
		sqlite3_finalize(pstmt);
		sqlite3_close_v2(db);
		exit(-1);
	}

	// fetch columns from our query
	while (sqlite3_step(pstmt) == SQLITE_ROW)
	{
		seed = string(reinterpret_cast<const char*>(sqlite3_column_text(pstmt, 0)));
		chromosome = string(reinterpret_cast<const char*>(sqlite3_column_text(pstmt, 1)));
		location = string(reinterpret_cast<const char*>(sqlite3_column_text(pstmt, 2)));
		three = string(reinterpret_cast<const char*>(sqlite3_column_text(pstmt, 3)));
		five = string(reinterpret_cast<const char*>(sqlite3_column_text(pstmt, 4)));
		score = string(reinterpret_cast<const char*>(sqlite3_column_text(pstmt, 5)));

		vector<string> chromosomeSplit = split(chromosome, delimeter);
		vector<string> locationSplit = split(location, delimeter);
		vector<string> threeSplit = split(three, delimeter);
		vector<string> fiveSplit = split(five, delimeter);
		vector<string> scoreSplit = split(score, delimeter);

		/* checks for if repeat is 3'/5'/both */
		if (threeSplit[0] == "" && fiveSplit[0] != "")
		{
			for (int i = 0; i < chromosomeSplit.size(); i++)
			{
				repeatSeqs += fiveSplit[i] + seed;
				repeatChroms.push_back(stoi(chromosomeSplit[i]));
				repeatLocations.push_back(stoll(locationSplit[i]));
				repeatScores.push_back(stoi(scoreSplit[i]));
			}
		}
		else if (threeSplit[0] != "" && fiveSplit[0] == "")
		{
			for (int i = 0; i < chromosomeSplit.size(); i++)
			{
				repeatSeqs += seed + threeSplit[i];
				repeatChroms.push_back(stoi(chromosomeSplit[i]));
				repeatLocations.push_back(stoll(locationSplit[i]));
				repeatScores.push_back(stoi(scoreSplit[i]));
			}
		}
		else if(threeSplit[0] != "" && fiveSplit[0] != "")
		{
			for (int i = 0; i < chromosomeSplit.size(); i++)
			{
				repeatSeqs += fiveSplit[i] + seed + threeSplit[i];
				repeatChroms.push_back(stoi(chromosomeSplit[i]));
				repeatLocations.push_back(stoll(locationSplit[i]));
				repeatScores.push_back(stoi(scoreSplit[i]));
			}
		}		
	}
	
	// 5 close the prepared statement
	sqlite3_finalize(pstmt);

	// 6 close the database
	sqlite3_close_v2(db);
}

/*
	function for parsing the input query sequences to score

	@param queryFilePath	=> file path to query file
	@param querySeqs		=> concatenated string of all sequences in the query file
	@param queryScores		=> vector of ints holding the scores of each query sequence
 */
void FileOperations::parseQueryFile(string &queryFilePath, string &querySeqs, vector<uint8_t> &queryScores)
{
	/* vars */
	string line;
	char queryFileDelimeter = ';';
	int queryCount = 0;
	fstream queryFile;
	
	/* open and verify file */
	queryFile.open(queryFilePath);
	if (queryFile.is_open())
	{
		while (getline(queryFile, line))
		{
			vector<string> lineSplit = split(line, queryFileDelimeter);
			querySeqs += lineSplit[1];
			queryScores.push_back(stoi(lineSplit[3]));
		}
		queryFile.close();
	}
	else
	{
		cerr << "Query input file was unable to be opened." << endl;
		exit(-1);
	}

	/* make sure query file had sequences */
	if (querySeqs.length() == 0)
	{
		cerr << "Query input file contained no sequences" << endl;
		exit(-1);
	}
}

/*
	function to open output file to write results to

	@param outputFilePath	=> file path for output file
	@param avgOutput		=> True: use average output, False: use detailed output format
 */
void FileOperations::openOutputFile(string &outputFilePath, bool &avgOutput)
{
	outputFile.open(outputFilePath);
	if (outputFile.is_open())
	{	
		/* Average output */
		if (avgOutput == true)
		{
			outputFile << "AVG OUTPUT" << endl;
		}
		/* Detailed output */
		else
		{
			outputFile << "DETAILED OUTPUT" << endl;
		}
	}
	else
	{
		cerr << "Output file couldn't be opened." << endl;
		exit(-1);
	}
}

/*

*/
void FileOperations::writeResults(bool &avgOutput, string &querySeq, vector<double> &uniqueScores, vector<unsigned long> &uniqueIndexes, vector<double> &repeatScores, vector<unsigned long> &repeatIndexes, vector<long long> &repeatLocations, vector<int> &repeatChroms, string &repeatSeqs, vector<long long> &uniqueLocations, vector<int> &uniqueChroms, string &uniqueSeqs, int &seqLength)
{
	double averageScore = accumulate(uniqueScores.begin(), uniqueScores.end(), 0.0) + accumulate(repeatScores.begin(), repeatScores.end(), 0.0);
	if (uniqueScores.size() + repeatScores.size() != 0)
	{
		averageScore /= (uniqueScores.size() + repeatScores.size());
	}

	outputFile << fixed;
	outputFile << setprecision(6);
	outputFile << querySeq << ":" << averageScore << endl;

	if (avgOutput == false)
	{
		for (int i = 0; i < uniqueScores.size(); i++)
		{
			outputFile <<uniqueScores[i] << "," << uniqueChroms[uniqueIndexes[i]] << "," << uniqueLocations[uniqueIndexes[i]] << "," << uniqueSeqs.substr(uniqueIndexes[i] * seqLength, seqLength) << endl;
		}
		for (int i = 0; i < repeatScores.size(); i++)
		{
			//cout << repeatSeqs.substr(repeatIndexes[i] * seqLength, seqLength) << endl;
			outputFile << repeatScores[i] << "," << repeatChroms[repeatIndexes[i]] << "," << repeatLocations[repeatIndexes[i]] << "," << repeatSeqs.substr(repeatIndexes[i] * seqLength, seqLength) << endl;
		}
	}
}

/*
	function to close output file
*/
void FileOperations::closeOutputFile()
{
	outputFile.close();
}

/*
	function split a given string based on a delimter

	@param s		=> string to split
	@param delimter	=> delimters to use for splitting

	@return splitStr	=> vector of substrings from split
 */
vector<string> FileOperations::split(string s, char &delimiter)
{
	size_t pos = 0;
	string token;
	vector <string> splitStr;
	while ((pos = s.find(delimiter)) != string::npos)
	{
		token = s.substr(0, pos);
		splitStr.push_back(token);
		s.erase(0, pos + 1);
	}
	splitStr.push_back(s);
	return splitStr;
}

/*
	function to load the default HSU MATRIX-spCas9-2013

	@param hsuMatrix	=> map to hold the HSU matrix data for OffTarget class
*/
void FileOperations::loadDefaultHsuMatrix(map<string, vector<double>> &hsuMatrix)
{
	hsuMatrix["GT"] = { 1, 1.6533,0.9030,1.5977,0.9235,0.8070,0.9632,1.0163,0.2658,0.7119,1.5211,0.6942,1.0434,0.5255,0.8981,0.7164,0.8399,0.5376,0.2821,0.6898 };
	hsuMatrix["AC"] = { 1, 1.5142,1.1597,1.6582,0.9924,0.0247,0.5522,1.8687,0.7737,0.9270,0.7292,0.4842,0.4824,0.7060,1.0221,0.0181,0.3496,0.1811,0.1362,0.2700 };
	hsuMatrix["GG"] = { 1, 1.3234,1.4157,1.2967,1.2060,0.9524,0.2304,1.0163,0.8100,1.1559,0.7075,1.5791,0.3490,0.0899,0.0497,0.0045,0.2267,0.2153,0.5250,0.4965 };
	hsuMatrix["TG"] = { 1, 1.5366,1.2771,1.2689,1.2197,1.0645,0.7791,1.2445,0.9885,1.5319,0.7184,1.7381,0.4166,0.1285,0.0720,0.0549,0.2261,0.3119,0.1343,0.0601 };
	hsuMatrix["TT"] = { 1, 1.7347,0.8215,1.1579,1.1816,0.7380,0.9004,0.8368,0.2997,0.6210,1.1400,0.3561,0.6192,0.1799,0.2665,0.2793,0.2613,0.1152,0.1680,0.3372 };
	hsuMatrix["CA"] = { 1, 1.0186,0.9649,2.2504,1.6222,0.2405,0.7561,1.0651,0.1102,1.4293,0.3533,0.6178,0.7269,0.1165,0.0367,0.5013,0.4147,0.1786,0.5315,0.1664 };
	hsuMatrix["CT"] = { 1, 1.6719,0.9688,1.0732,1.0869,0.6475,1.0142,0.8635,0.3059,0.4487,0.9046,0.4327,0.5576,0.1379,0.0722,0.3279,0.2420,0.0433,0.1351,0.4403 };
	hsuMatrix["GA"] = { 1, 1.1662,0.4544,2.7867,1.0461,0.6036,0.8132,0.7875,0.6882,1.3655,0.1240,0.1953,0.2497,0.0132,0.0227,0.0478,0.3682,0.3175,0.5621,0.4588 };
	hsuMatrix["AA"] = { 1, 1.1916,1.0954,2.8709,1.1310,0.5160,0.6439,1.0322,0.5356,1.2868,0.0780,0.2592,0.2675,0.0469,0.0252,0.0052,0.0218,0.1718,0.6970,0.2720 };
	hsuMatrix["AG"] = { 1, 1.4786,1.0820,1.2952,0.7450,0.9763,0.4912,0.9272,0.6022,1.0375,0.3047,0.8210,0.0303,0.0365,0.0592,0.0253,0.1553,0.1006,0.2175,0.0275 };
	hsuMatrix["TC"] = { 1, 1.0400,0.9954,1.6466,1.3410,0.0102,0.5428,2.3401,0.4367,0.2143,0.3405,0.2640,0.0935,0.0462,0.0688,0.0165,0.3659,0.0546,0.0857,0.2543 };
	hsuMatrix["CC"] = { 1, 1.0345,1.0478,1.0507,1.4075,0.0540,0.6396,2.0810,0.4585,0.1555,0.1369,0.1026,0.0417,0.0105,0.0458,0.0099,0.2114,0.0552,0.0253,0.0596 };
}