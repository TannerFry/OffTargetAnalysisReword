#include "OffTarget.h"

#include <chrono>
#include <fstream>
#include "Score.h"

using namespace std;
using std::chrono::high_resolution_clock;
using std::chrono::duration_cast;
using std::chrono::duration;
using std::chrono::milliseconds;
using std::chrono::seconds;

int main(int argc, char *argv[])
{	
	auto start = high_resolution_clock::now();

	/* create OffTarget object */
	OffTarget OT;

	/* store input arguments */
	cout << "Parsing Input Arguments" << endl;
	OT.parseInputArguments(argc, argv);

	/* parse data needed for algorithm */
	cout << "Loading data for algorithm" << endl;
	OT.getAlgorithmData();

	/* run the OffTarget algorithm */
	cout << "Running OffTarget Analysis" << endl;
	OT.run();

	auto end = high_resolution_clock::now();
	
	cout << "Elapsed time in seconds: "
		<< chrono::duration_cast<chrono::seconds>(end - start).count()
		<< " sec" << endl;

	//system("pause");
	return 0;
}

