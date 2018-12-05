/**************************************************************************************************************************************
Description: 
        Combine the carbon nanotube and carbon nanosheet

Date:
         12/05/2018
		 
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

string Cnt_Coor_File="cnt.xyz";
string Sheet_Coor_File="sheet.xyz";
string Out_File="cnt_sheet.xyz";
float Kvalue=2.0;

float KKvalue=1.4;

int printhelp(){
	cout << "\t  -i CNT coordinate file name" << endl;
	cout << "\t  -f CNS coordinate file name" << endl;
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
        case 'i': Cnt_Coor_File = argv[i+1];break;
		case 'f': Sheet_Coor_File = argv[i+1];break;
		case 'o': Out_File = argv[i+1];break;
		case 'b': Kvalue = atof(argv[i+1]);break;
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

vector< vector<float> >TransBy_Center(vector< vector<float> > Coor){

    float center_x=0.0;
	float center_y=0.0;
	float center_z=0.0;
	
	for(int i=0; i<Coor.size(); i++){
		center_x+=Coor[i][0];
		center_y+=Coor[i][1];
		center_z+=Coor[i][2];
	}
	vector< vector<float> > Center_Coor;
    center_x=center_x/Coor.size();
	center_y=center_y/Coor.size();
	center_z=center_z/Coor.size();
	
	for(int i=0; i<Coor.size(); i++){
		vector<float> Per_Atom;
		Per_Atom.push_back(Coor[i][0]-center_x);
		Per_Atom.push_back(Coor[i][1]-center_y);
		Per_Atom.push_back(Coor[i][2]-center_z);
		Center_Coor.push_back(Per_Atom);
	}
	return Center_Coor;
}

 int main(int argc, char * argv[]){
	 Parse_Para(argc, argv);
	 vector<string> Atom_CNT, Atom_CNS;
	 vector< vector<float> > Coor_CNT, Coor_CNS;
	 ofstream outfile(Out_File.c_str(), ofstream::out);
/* Read in the XYZ files, get the number of atoms and types*/
	vector< vector<string> > Inform_XYZ;
	Inform_XYZ=Load_XYZ_File(Cnt_Coor_File.c_str());
	int natoms;
	natoms= atoi(Inform_XYZ[0][0].c_str());
	for(int i=0; i<natoms; i++){
		int k=i+1;
		vector<float> Per_Vector;
		for(int j=0; j<Inform_XYZ[k].size(); j++){
			if(j==0){
				Atom_CNT.push_back(Inform_XYZ[k][j]);
			}else{
				Per_Vector.push_back(atof(Inform_XYZ[k][j].c_str()));
			}
		}
		Coor_CNT.push_back(Per_Vector);
	}
	
	Inform_XYZ.clear();
	Inform_XYZ=Load_XYZ_File(Sheet_Coor_File.c_str());
	natoms = atoi(Inform_XYZ[0][0].c_str());
	for(int i=0; i<natoms; i++){
		int k=i+1;
		vector<float> Per_Vector;
		for(int j=0; j<Inform_XYZ[k].size(); j++){
			if(j==0){
				Atom_CNS.push_back(Inform_XYZ[k][j]);
			}else{
				Per_Vector.push_back(atof(Inform_XYZ[k][j].c_str()));
			}
		}
		Coor_CNS.push_back(Per_Vector);
	}
	
	Coor_CNT=TransBy_Center(Coor_CNT);
	Coor_CNS=TransBy_Center(Coor_CNS);
	
	vector<float>  Min_CNT, Max_CNT;
	for(int i=0; i<Coor_CNT.size(); i++){
			for(int j=0; j<Coor_CNT[i].size(); j++){
				if(i==0){
				Min_CNT.push_back(Coor_CNT[i][j]);
				Max_CNT.push_back(Coor_CNT[i][j]);
				}else{
				float temp_min=(Min_CNT[j]<Coor_CNT[i][j])?Min_CNT[j]:Coor_CNT[i][j];
				Min_CNT[j]=temp_min;
				float temp_max=(Max_CNT[j]>Coor_CNT[i][j])?Max_CNT[j]:Coor_CNT[i][j];
				Max_CNT[j]=temp_max;
				}
			}
	}
	
	if(Max_CNT[0]!=Max_CNT[1]){
		cerr<<"The nanotube should be distributed along Z axis"<<endl;
	}

	vector<int> CNS_Print;
	
	for(int i=0; i<Coor_CNS.size(); i++){
		CNS_Print.push_back(1);
		for(int k=0; k<Coor_CNT.size(); k++){
			float dx=Coor_CNS[i][0]-Coor_CNT[k][0];
			float dy=Coor_CNS[i][1]-Coor_CNT[k][1];
			float dz=Coor_CNS[i][2]-Coor_CNT[k][2];
			float dist=sqrt(dx*dx+dy*dy+dz*dz);
			if(dist<KKvalue){
				CNS_Print[i]=0;
			}
		}
		
		if(abs(Coor_CNS[i][0])<Max_CNT[0]&&abs(Coor_CNS[i][1])<Max_CNT[1]){
			CNS_Print[i]=0;
		}
	}
	
	for(int i=0; i<Coor_CNS.size(); i++){
		if(CNS_Print[i]==1){
			outfile<<Atom_CNS[i]<<"\t"<<Coor_CNS[i][0]<<"\t"<<Coor_CNS[i][1];
			outfile<<"\t"<<Min_CNT[2]+Kvalue<<endl;
		}
	}
	
	for(int i=0; i<Coor_CNS.size(); i++){
		if(CNS_Print[i]==1){
			outfile<<Atom_CNS[i]<<"\t"<<Coor_CNS[i][0]<<"\t"<<Coor_CNS[i][1];
			outfile<<"\t"<<Max_CNT[2]-Kvalue<<endl;
		}
	}
	
	for(int i=0; i<Coor_CNT.size(); i++){
		outfile<<Atom_CNT[i];
		for(int k=0; k<Coor_CNT[i].size(); k++){
			outfile<<"\t"<<Coor_CNT[i][k];
		}
		outfile<<endl;
	}
	
 }
