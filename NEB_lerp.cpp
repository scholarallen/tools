/**************************************************************************************************************************************
Description: 
        Predict coordinates used for the calculation of Nudged Elastic Band(NEB) by Linear 
Interpolartion Algorithm. This code only can be used for isolated molecules but can not be 
used for the crystal system.

Date:
         11/27/2018
		 
Authors:
         Allen Wang 

Any problems, please feel free to let me know. Thanks!
E-mail: 
         scholar_allen@yeah.net
***************************************************************************************************************************************/

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <map>
#include <set>
#include <ctime>
#include <math.h>
#include <string.h>
#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <iomanip>

using namespace std;

string Initial_Coor_File="IC.xyz";
string Final_Coor_File="FC.xyz";
string Out_File="Output_Coor_";
int binvalue=1;

int printhelp(){
	cout << "\t  -i Reaction coordinate file name" << endl;
	cout << "\t  -f Production coordinate file name" << endl;
	cout << "\t  -o Output file name" << endl;
	cout << "\t  -b Bins" << endl;
	cout << "\t  -h Help" << endl;
}

int Parse_Para(int argc, char * argv[]){
	
	if (argc ==1) 
		printhelp();
	
	int i=1;
	while(i<argc){
		if(argv[i][0]!='-'){
			cerr << "Argument # " << i << " Error : Arguments must start with -" << endl;
			exit(0);
		}
		switch(argv[i][1]){
        case 'i': Initial_Coor_File = argv[i+1];break;
		case 'f': Final_Coor_File = argv[i+1];break;
		case 'o': Out_File = argv[i+1];break;
		case 'b': binvalue = atoi(argv[i+1]);break;
		case 'h': printhelp(); break;
		}
		i=i+2;
	}
	
	return 0;
}

vector< vector<string> > Load_XYZ_File(const char * infilename){
	ifstream infile(infilename, ifstream::in);
    if(!infile){
                 cerr << "Error: Cannot open file : " << infilename << endl;
                 exit(0);
                 }
	vector< vector<string> > return_vec;
	string buffer, temp;
	
	getline(infile, buffer);
	vector<string> temp_vec;
	temp_vec.push_back(buffer);
	return_vec.push_back(temp_vec);

	getline(infile, buffer);
	while(getline(infile, buffer)){
		stringstream strin(buffer);
		temp_vec.clear();
		while(strin>>temp){
			temp_vec.push_back(temp);
		}
		return_vec.push_back(temp_vec);
	}
	
	return return_vec;
}
 
 
string intTstring(int input_int){
	ostringstream os;
	os<<input_int;
	return os.str();
}

void Print_Results(vector< vector<string> > ICoor,vector< vector<float>  >Diffs, int nframe){
	string outname=Out_File+intTstring(nframe)+".xyz";
	ofstream outfile(outname.c_str(), ofstream::out);
	outfile<<atoi(ICoor[0][0].c_str())<<endl;
	outfile<<endl;
	for(int i=0; i<Diffs.size(); i++){
		int k=i+1;
		outfile<<setw(2)<<ICoor[k][0]<<setw(17)<<setprecision(8)<<atof(ICoor[k][1].c_str())+Diffs[i][0]*i/binvalue;
		outfile<<setw(17)<<setprecision(8)<<atof(ICoor[k][2].c_str())+Diffs[i][1]*i/binvalue;
		outfile<<setw(17)<<setprecision(8)<<atof(ICoor[k][3].c_str())+Diffs[i][2]*i/binvalue<<endl;
	}
	outfile<<endl;
}
 
 int main(int argc, char * argv[]){
	 Parse_Para(argc, argv);
/* Read in the XYZ files and make sure that the number of atoms and types are the same*/
     vector< vector<string> > Contex_I, Contex_F;
	 Contex_I=Load_XYZ_File(Initial_Coor_File.c_str());
	 Contex_F=Load_XYZ_File(Final_Coor_File.c_str());
	 
	 if(atoi(Contex_I[0][0].c_str())!=atoi(Contex_F[0][0].c_str()) || Contex_I.size()!=Contex_F.size()){
		 cerr<<"WARNING:: Atom numbers are not the same, hope you know what you are doing"<<endl;
	 }
	 
	 vector< vector<float> > Max_Diff;
	 for(int i=0; i<atoi(Contex_I[0][0].c_str()); i++){
		 int k=i+1;
		 vector<float> Diff_Temp;
		 Diff_Temp.push_back(atof(Contex_F[k][1].c_str())-atof(Contex_I[k][1].c_str()));
		 Diff_Temp.push_back(atof(Contex_F[k][2].c_str())-atof(Contex_I[k][2].c_str()));
		 Diff_Temp.push_back(atof(Contex_F[k][3].c_str())-atof(Contex_I[k][3].c_str()));
		 Max_Diff.push_back(Diff_Temp);
	 }
	 
	 for(int i=0; i<binvalue+1; i++){
		 Print_Results(Contex_I, Max_Diff, i);
	 }
 }
 
 /* */
 