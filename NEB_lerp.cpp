/**************************************************************************************************************************************
Description: 
        Predict coordinates used for the calculation of Nudged Elastic Band(NEB) by Linear 
Interpolartion Algorithm. This code only can be used for isolated molecules but can not be 
used for the crystal system.

Date:
         11/27/2018
		 05/16/2019 *.gjf format
		 
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
#include <string>
#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <iomanip>

using namespace std;

string Test_Name;
string Initial_Coor_File="IC.xyz";
string Final_Coor_File="FC.xyz";
string Out_File="Output_Coor_";
vector<string> GJFA, GJFB;
int Start_Sign=0;
int binvalue=1;

int printhelp(){
	cout << "\t  -i Reaction coordinate file name" << endl;
	cout << "\t  -f Production coordinate file name" << endl;
	cout << "\t  -s Start Serial" << endl;
	cout << "\t  -o Output file name" << endl;
	cout << "\t  -b Bins" << endl;
	cout << "\t  -h Help" << endl;
	exit(0);
	return 0;
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
        case 'i': Initial_Coor_File = argv[i+1];Test_Name=argv[i+1];break;
		case 'f': Final_Coor_File = argv[i+1];break;
		case 's': Start_Sign = atoi(argv[i+1]);break;
		case 'o': Out_File = argv[i+1];break;
		case 'b': binvalue = atoi(argv[i+1]);break;
		case 'h': printhelp(); break;
		}
		i=i+2;
	}
	
	return 0;
}

vector<string> Filter_string(string Input_String, const char * split_char){
	
	vector<string> real_string;
	
	char *pr=strtok((char *)Input_String.c_str(),split_char);
	while(pr!=NULL){
		real_string.push_back(pr);
		pr=strtok(NULL,split_char);
	}
	return real_string;
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

void trim(string &s){

    int index = 0;
    if( !s.empty()){
        while( (index = s.find(' ',index)) != string::npos){
            s.erase(index,1);
        }
    }

}

vector< vector<string> > Load_GJF_File(const char * infilename){
	ifstream infile(infilename, ifstream::in);
    if(!infile){
                 cerr << "Error: Cannot open file : " << infilename << endl;
                 exit(0);
                 }
	vector< vector<string> > return_vec;
	string buffer, temp;
	int Nempty=0;
	bool Read_Coor=false;
        bool skipR=true;
	GJFA.clear();
	GJFB.clear();
	while(getline(infile, buffer)){
		string buffer1=buffer;
		trim(buffer1);
		if(buffer1.length()==0){
			Nempty++;
		}
		if(Nempty<2){
			GJFA.push_back(buffer);
		}
		if(Nempty==2&&skipR){
			GJFA.push_back(buffer);
			getline(infile, buffer);
			GJFA.push_back(buffer);
			getline(infile, buffer);
			//GJFA.push_back(buffer);
			Read_Coor=true;
                        skipR=false;
		}
		if(Nempty>=3){
			Read_Coor=false;
			GJFB.push_back(buffer);
		}
		if(Read_Coor){
			stringstream strin(buffer);
			vector<string> temp_vec;
			while(strin>>temp){
				temp_vec.push_back(temp);
			}
			return_vec.push_back(temp_vec);
		}
		
	}

	return return_vec;
}

string intTstring(int input_int){
	ostringstream os;
	os<<input_int;
	return os.str();
}

void Print_Results(vector< vector<string> > ICoor,vector< vector<float>  >Diffs, int nframe,string Ptype){
	string outname;
	int Isign=Start_Sign+nframe;
	if(Ptype=="xyz"){
		outname=Out_File+intTstring(Isign)+".xyz";
	}
	if(Ptype=="gjf"){
		outname=Out_File+intTstring(Isign)+".gjf";
	}
	ofstream outfile(outname.c_str(), ofstream::out);
	if(Ptype=="xyz"){
		outfile<<atoi(ICoor[0][0].c_str())<<endl;
		outfile<<endl;
	}
	if(Ptype=="gjf"){
		for(int i=0; i<GJFA.size(); i++){
			outfile<<GJFA[i]<<endl;
		}
	}
	for(int i=0; i<Diffs.size(); i++){
		 int k;
		 if(Ptype=="xyz"){
			k=i+1;
		 }
		 if(Ptype=="gjf"||Ptype=="com"){
			k=i;
		 }
		outfile<<setw(2)<<ICoor[k][0]<<setw(17)<<setprecision(8)<<atof(ICoor[k][1].c_str())+Diffs[i][0]*nframe/binvalue;
		outfile<<setw(17)<<setprecision(8)<<atof(ICoor[k][2].c_str())+Diffs[i][1]*nframe/binvalue;
		outfile<<setw(17)<<setprecision(8)<<atof(ICoor[k][3].c_str())+Diffs[i][2]*nframe/binvalue<<endl;
	}
	if(Ptype=="gjf"){
		for(int i=0; i<GJFB.size(); i++){
			outfile<<GJFB[i]<<endl;
		}
	}
	outfile<<endl;
}
 
 int main(int argc, char * argv[]){
	 Parse_Para(argc, argv);
	 vector< vector<string> > Contex_I, Contex_F;
	 vector<string> File_Type_vec=Filter_string(Test_Name,".");
	 int vec_size=File_Type_vec.size()-1;
	 int natoms;
	 if(File_Type_vec[vec_size]=="xyz"){
/* Read in the XYZ files and make sure that the number of atoms and types are the same*/
		Contex_I=Load_XYZ_File(Initial_Coor_File.c_str());
		Contex_F=Load_XYZ_File(Final_Coor_File.c_str());
		natoms=atoi(Contex_I[0][0].c_str());
	 }
	 if(File_Type_vec[vec_size]=="gjf"||File_Type_vec[vec_size]=="com"){
		Contex_I=Load_GJF_File(Initial_Coor_File.c_str());
		Contex_F=Load_GJF_File(Final_Coor_File.c_str());
		natoms=Contex_I.size();
	 }

	 
	 if(atoi(Contex_I[0][0].c_str())!=atoi(Contex_F[0][0].c_str()) || Contex_I.size()!=Contex_F.size()){
		 cerr<<"WARNING:: Atom numbers are not the same, hope you know what you are doing"<<endl;
	 }
	 vector< vector<float> > Max_Diff;
	 for(int i=0; i<natoms; i++){
		 int k;
		 if(File_Type_vec[vec_size]=="xyz"){
			k=i+1;
		 }
		 if(File_Type_vec[vec_size]=="gjf"||File_Type_vec[vec_size]=="com"){
			k=i;
		 }
		 vector<float> Diff_Temp;
		 Diff_Temp.push_back(atof(Contex_F[k][1].c_str())-atof(Contex_I[k][1].c_str()));
		 Diff_Temp.push_back(atof(Contex_F[k][2].c_str())-atof(Contex_I[k][2].c_str()));
		 Diff_Temp.push_back(atof(Contex_F[k][3].c_str())-atof(Contex_I[k][3].c_str()));
		 Max_Diff.push_back(Diff_Temp);
	 }
	 
	 for(int i=0; i<binvalue+1; i++){
		 
		if(File_Type_vec[vec_size]=="xyz"){
			Print_Results(Contex_I, Max_Diff, i, "xyz");
		}
		if(File_Type_vec[vec_size]=="gjf"||File_Type_vec[vec_size]=="com"){
			Print_Results(Contex_I, Max_Diff, i, "gjf");
		}
	 }
}
 
 /* */
 
