# **CASPER Off-Target Algorithm**

This algorithm uses CASPER's indexes from Seq-Finder to find and score potential off-target sites for input gRNAs.

## Download And Compile Sqlite3
### Linux: 
1. Open terminal
2. Run the following command: `sudo apt-get install sqlite3`

### Mac and Linux (if you did not use apt-get):
1. Download sqlite3 source code
2. Open terminal
3. CD to sqlite3 source code folder
3. Build using the following command: `gcc -c sqlite3.c`
	* The build command should generate a .o file that will be used with compiling OT.
4. Copy the .o file to the OT source code folder

### Windows (Visual Studio 2017):
1. Download sqlite3 source code
2. Open "Developer Command Prompt for VS 2017"
3. CD to sqlite3 source code folder
4. Build the library files by running the following command: `lib /DEF:sqlite3.def /OUT:sqlite3.lib /MACHINE:x64`

## Download and Compile OT
### Linux (if you used apt-get install sqlite3):
1. Download OT source code for Linux (in Repository)
2. Open terminal
3. CD to OT source code folder
3. Run Command to compile OT: `g++ -std=c++11 *.cpp -pthread -lsqlite3 -o OT`

### Mac and Linux (if you manually built sqlite3 .o file, make sure sqlite3 .o file is in same folder as OT source code):
1. Download OT source code for Mac or Linux (in Repository)
2. Open terminal
3. CD to OT source code folder
3. Run command to compile OT: `g++ -std=c++11 *.cpp -pthread sqlite3.o -o OT`

### Windows (Visual Studio 2017):
1. Download OT source for Windows (in Repository)
2. Open Visual Studio 2017
3. Create a New Project
4. Import the OT source code files
5. Go to Project->Properties and do the following:
	* Set Configuration to "All Configurations"
	* Set Platform to "All Platforms"
	* In Configuration Properties->VC++ Directories add `C:\Path\To\sqlite3;` to "Include Directories"
	* In Configuration Properties->VC++ Directories add `C:\Path\To\sqlite;` to "Include Libraries"
	* In C/C++->General add `C:\Path\To\sqlite3;` to "Additional Include Directories" and "Additional #using Directories"
	* In Linker->General add `C:\Path\To\sqlite;` to "Additional Library Directories"
	* In Linker->Input add `C:\Path\To\sqlite\sqlite3.lib;` in "Additonaly Dependencies
	* Make sure when you add these paths that there are ';' seperating all paths/object in the line.
6. For debugging, make sure Debug and x64 are selected before running
	* If debugging, make sure you set the debugging command arguments. See "How to run OT" below.
7. For compiling, make sure Release and x64 are selected before running
	

## How to run OT
* CD to the directory containing the OT executable
	* The command line arguments for OT are as follows: `query_file_path endonuclease cspr_file_path db_file_path output_file_path CASPERinfo_file_path max_num_mismatches threshold detailed_output_bool avg_output_bool hsu_matrix_name`

* Example command: `./OT query.txt asCas12 myfile_asCas12.cspr myfile_asCas12_repeats.db output.txt CASPERinfo 5 0.05 TRUE FALSE "MATRIX:HSU MATRIX-asCas12-2016"`
