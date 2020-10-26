/**************************************************************************************************************************************
Description: 
        Calculate the IP-correlation function of one type vector(formed by two molecules or atom) change with time. 
These molcules (atom) and its mass should be given by an index file. 
eg.
###Pair #1
C1  H1  O1
12.0 1.0 16.0
###Pair #2
N1  H2
14.0 1.0
As shown in above table, the code will calculate the vector that were formed by Pair #1 and Pair #2. 

##Usage::
    g++  tools_ilp.cpp include/libxdrfile.a -I include/

Date:
         09/24/2020
		 
Authors:
         Allen Wang, Guokui Liu
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
#include <omp.h>
#include <math.h>
#include <unistd.h>
#include <stdlib.h>

extern "C" {
#include "include/xdrfile.h"
#include "include/xdrfile_trr.h"
#include "include/xdrfile_xtc.h"
}

using namespace std;

int natoms=0;
string coord_name="traj.xtc";
string xyzfile_name="conf.xyz";
string indexfile_name="index.ndx";
int Max_Corr_Time=100;
vector<string> Inform_Pair1, Inform_Pair2;
map<int,string> Inform_Atom_Name;
map<string, float> Inform_Mass;

int printhelp(){
    cout << "\t calculate the IP-correlation function" << endl;
    cout << "Usage: " << endl;
    cout << "option Value" << endl;
    cout << "Options: " << endl;
    cout << "\t  -f trajectory file name (.xtc)"<<endl;
	cout << "\t  -p reference file name (.xyz)"<<endl;
    cout << "\t  -n index file name (.ndx)"<<endl;
    cout << "\t  -t Max correlation time, default 100ps"<<endl;
    cout << "\t  -h Help" << endl;
    exit(0);
    return 0;
    }

void Parse_Para(int argc, char * argv[]){

	int i=1;
	while(i<argc){
		if(argv[i][0]!='-'){
			cerr << "Argument # " << i << " Error : Arguments must start with -" << endl;
			exit(0);
		};
		switch(argv[i][1]){
        case 'f': coord_name = argv[i+1];break;
        case 'n': indexfile_name = argv[i+1];break;
        case 'p': xyzfile_name = argv[i+1];break;
        case 't': Max_Corr_Time = atoi(argv[i+1]);break;
		case 'h': printhelp();
		default : cerr << "Error: Unrec argument " << argv[i] << endl;
                  printhelp();
			      break;
		}
	 	i=i+2;
	}
	if (indexfile_name.size() == 0){
        cerr << printhelp();
        exit(0);
    }else{
		ifstream infile(indexfile_name.c_str(), ifstream::in);
		if(!infile){
			cerr << printhelp();
			exit(0);
		}else{
			infile.close();
		}
	}
}

int GetIndex(const char* indexfile){
	
	ifstream infile(indexfile, ifstream::in);
    if(!infile){
        cout<<"Error: Cannot open file : " << indexfile << endl;
		return 0;
   }
	
	string buffer, temp;
	
	getline(infile,buffer);
	stringstream strin1(buffer);
	while(strin1>>temp){
		Inform_Pair1.push_back(temp);
	}
	getline(infile,buffer);
	stringstream strin1m(buffer);
	int m1=0;
	while(strin1m>>temp){
		Inform_Mass[Inform_Pair1[m1]]=atof(temp.c_str());
		m1++;
	}
	getline(infile,buffer);
	stringstream strin2(buffer);
	while(strin2>>temp){
		Inform_Pair2.push_back(temp);
	}
	getline(infile,buffer);
	stringstream strin2m(buffer);
	int m2=0;
	while(strin2m>>temp){
		Inform_Mass[Inform_Pair2[m2]]=atof(temp.c_str());
		m2++;
	}
	return 0;
}

int GetXyz(const char* xyzfile){
	ifstream infile(xyzfile, ifstream::in);
    if(!infile){
        cout<<"Error: Cannot open file : " << xyzfile << endl;
		return 0;
   }
	string buffer, temp;
	
	getline(infile,buffer);
	stringstream strin1(buffer);
	strin1>>temp;
	getline(infile,buffer);
	natoms=atoi(temp.c_str());
	for(int i=0; i<natoms; i++){
		getline(infile,buffer);
		stringstream strin(buffer);
		strin>>temp;
		Inform_Atom_Name[i]=temp;
	}
	
	return 0;
}

int main(int argc, char *argv[]){
	Parse_Para(argc, argv);
	GetXyz(xyzfile_name.c_str());
	GetIndex(indexfile_name.c_str());
	rvec *x;
    matrix box;
    XDRFILE *xtc;
	
    int step;
    float time,p;
    char * filename;
    filename=(char*) coord_name.c_str();
    vector<float> Inform_Time;
	
    xtc=xdrfile_open(filename,"r");
    int read_return=read_xtc_natoms(filename,&natoms);
    x=(rvec * )calloc(natoms,sizeof(x[0]));
	
	map<string,string> Pair1_Cont, Pair2_Cont;
	for(int i=0; i<Inform_Pair1.size(); i++){
		Pair1_Cont[Inform_Pair1[i]]="T";
	}
	
	for(int i=0; i<Inform_Pair2.size(); i++){
		Pair2_Cont[Inform_Pair2[i]]="T";
	}
	
	vector< vector<float> > Vecx_Pair1, Vecy_Pair1, Vecz_Pair1;
	vector< vector<float> > Vecx_Pair2, Vecy_Pair2, Vecz_Pair2;
	while(true){
		read_return=read_xtc(xtc,natoms,&step,&time,box,x,&p);
        if(read_return!=0) break;
		Inform_Time.push_back(time);
		vector<float> PerPair1_x, PerPair1_y, PerPair1_z;
		vector<float> PerPair2_x, PerPair2_y, PerPair2_z;
		for(int i=0; i<natoms; i++){
			if(Pair1_Cont[Inform_Atom_Name[i]]=="T"){
				PerPair1_x.push_back(x[i][0]);
				PerPair1_y.push_back(x[i][1]);
				PerPair1_z.push_back(x[i][2]);
			}
			if(Pair2_Cont[Inform_Atom_Name[i]]=="T"){
				PerPair2_x.push_back(x[i][0]);
				PerPair2_y.push_back(x[i][1]);
				PerPair2_z.push_back(x[i][2]);
			}
		}
		int npairs1=PerPair1_x.size()/Inform_Pair1.size();
		int nt1=0;
		vector<float> Pair1_MassCenter_x, Pair1_MassCenter_y, Pair1_MassCenter_z;
		for(int i=0; i<npairs1; i++){
			float txt1x=0;
			float txt1y=0;
			float txt1z=0;
			float tmt1=0;
			for(int k=0; k<Inform_Pair1.size(); k++){
				txt1x+=Inform_Mass[Inform_Pair1[k]]*PerPair1_x[nt1];
				txt1y+=Inform_Mass[Inform_Pair1[k]]*PerPair1_y[nt1];
				txt1z+=Inform_Mass[Inform_Pair1[k]]*PerPair1_z[nt1];
				tmt1+=Inform_Mass[Inform_Pair1[k]];
				nt1++;
			}
			Pair1_MassCenter_x.push_back(txt1x/tmt1);
			Pair1_MassCenter_y.push_back(txt1y/tmt1);
			Pair1_MassCenter_z.push_back(txt1z/tmt1);
		}
		
		int npairs2=PerPair2_x.size()/Inform_Pair2.size();
		nt1=0;
		vector<float> Pair2_MassCenter_x, Pair2_MassCenter_y, Pair2_MassCenter_z;
		for(int i=0; i<npairs2; i++){
			float txt2x=0;
			float txt2y=0;
			float txt2z=0;
			float tmt2=0;
			for(int k=0; k<Inform_Pair2.size(); k++){
				txt2x+=Inform_Mass[Inform_Pair2[k]]*PerPair2_x[nt1];
				txt2y+=Inform_Mass[Inform_Pair2[k]]*PerPair2_y[nt1];
				txt2z+=Inform_Mass[Inform_Pair2[k]]*PerPair2_z[nt1];
				tmt2+=Inform_Mass[Inform_Pair2[k]];
				nt1++;
			}
			Pair2_MassCenter_x.push_back(txt2x/tmt2);
			Pair2_MassCenter_y.push_back(txt2y/tmt2);
			Pair2_MassCenter_z.push_back(txt2z/tmt2);
		}
		
		Vecx_Pair1.push_back(Pair1_MassCenter_x);
		Vecx_Pair2.push_back(Pair2_MassCenter_x);
		Vecy_Pair1.push_back(Pair1_MassCenter_y);
		Vecy_Pair2.push_back(Pair2_MassCenter_y);
		Vecz_Pair1.push_back(Pair1_MassCenter_z);
		Vecz_Pair2.push_back(Pair2_MassCenter_z);
	}
	
	vector< vector<int> > Num_Matrix;
	for(int i=0; i<Vecx_Pair1.size(); i++){
		vector<int> Per_Matrix;
		for(int k=0; k<Vecx_Pair1[i].size(); k++){
			float mindist=10000.0;
			int n_mindist=0;
			for(int g=0; g<Vecx_Pair2[i].size(); g++){
				float dx=Vecx_Pair1[i][k]-Vecx_Pair2[i][g];
				float dy=Vecy_Pair1[i][k]-Vecy_Pair2[i][g];
				float dz=Vecz_Pair1[i][k]-Vecz_Pair2[i][g];
				float dist=sqrt(dx*dx+dy*dy+dz*dz);
				if(dist<mindist){
					mindist=dist;
					n_mindist=g;
				}
			}
			Per_Matrix.push_back(n_mindist);
		}
		Num_Matrix.push_back(Per_Matrix);
	}
	float Delt_Time=Inform_Time[1]-Inform_Time[0];
	int Max_Corr_Step=ceil(Max_Corr_Time/Delt_Time);
	
	for(int i=1; i<Max_Corr_Step+1; i++){
		int nCycle=0;
		float nPartical=0;
		for(int k=0; k<Vecx_Pair1.size()-i; k++){
			for(int g=0; g<Vecx_Pair1[k].size(); g++){
				if(Num_Matrix[k][g]==Num_Matrix[k+i][g]){
					nCycle++;
				}
				nPartical++;
			}
		}
		if(nPartical!=0)
			cout<<i*Delt_Time<<"\t"<<nCycle/nPartical<<endl;
	}
	
}