/**************************************************************************************************************************************
Description: 
        Calculate the autocorrelation function of one type vector(formed by two molecules) change with time. These molcules should be 
given by an index file and its format similar with GROMACS format. 

eg.
1 3 5 7 9
2 4 6 8 10
As shown in above table, the code will calculate the vector that were formed by (1,2), (3,4), ...., (9,10). 

Date:
         05/01/2018
		 
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

string traj_file;
string index_file="index.ndx";
string output_file="ac.dat";
string first_frame="0";
string last_frame="2147483647";

vector< vector<float> > coor_data_vec1_x,coor_data_vec1_y,coor_data_vec1_z;
vector< vector<float> > coor_data_vec2_x,coor_data_vec2_y,coor_data_vec2_z;
vector< vector<int> > mol_indx;
float dt=0;
float time_last=0.0;
float time_first=0.0;
int Threadnum=0;
int nspin=1;

int printhelp(){
    cout << "\t calculate the auto-correlation function of angles" << endl;
    cout << "Usage: " << endl;
    cout << "option Value" << endl;
    cout << "Options: " << endl;
    cout << "\t  -f trajectory file name (.xtc)"<<endl;
    cout << "\t  -n index file name (.ndx)"<<endl;
    cout << "\t  -b First frame (ps) to read from trajectory"<<endl;
    cout << "\t  -e Last frame (ps) to read from trajectory"<<endl;
    cout << "\t  -t Cpu core number, default is auto" <<endl;
    cout << "\t  -p print bin, default is 1" <<endl;
    cout << "\t  -o output file name, default is 'ac.dat' "<<endl;
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
        case 'f': traj_file = argv[i+1];break;
        case 'n': index_file = argv[i+1];break;
        case 'b': first_frame = argv[i+1];break;
        case 'e': last_frame = argv[i+1];break;
        case 't': Threadnum = atoi(argv[i+1]);break;
        case 'p': nspin = atoi(argv[i+1]);break;
        case 'o': output_file = argv[i+1];break;
		case 'h': printhelp();
		default : cerr << "Error: Unrec argument " << argv[i] << endl;
                  printhelp();
			      break;
		}
		i=i+2;
	}
	if (traj_file.size() == 0){ 
                   cerr << printhelp();
                   exit(0);
                   }
}


void load_index(const char* infilename){
     ifstream infile(infilename, ifstream::in);
     if(!infile){
                 cout<<"Error: Cannot open file : " << infilename << endl;
                 }
     string buffer;
     while(getline(infile,buffer)){
                                   stringstream strin(buffer);
                                   string temp;
                                   vector<int> per_indx;
                                   while(strin>>temp){
                                                      per_indx.push_back(atoi(temp.c_str()));
                                                      }
                                   mol_indx.push_back(per_indx);
                                   }
     
     }
     
void load_coor(string infilename){
     rvec *x;
     matrix box;
     XDRFILE *xtc;
    
     int natoms,step;
     float time,p;
     char * filename;
     filename=(char*) infilename.c_str();
    
     xtc=xdrfile_open(filename,"r");
     int read_return=read_xtc_natoms(filename,&natoms);
     x=(rvec * )calloc(natoms,sizeof(x[0]));
     
     int nframe=0;
     
     while(1){
              read_return=read_xtc(xtc,natoms,&step,&time,box,x,&p);
              if(read_return!=0) break;
              if(time>=atoi(first_frame.c_str()) && time<=atoi(last_frame.c_str())){
              if(nframe==0) time_first=time;
              vector<float> per_frame_x, per_frame_y, per_frame_z;
              for(int i=0; i<mol_indx[0].size(); i++){
                      per_frame_x.push_back(x[mol_indx[0][i-1]][0]);
                      per_frame_y.push_back(x[mol_indx[0][i-1]][1]);
                      per_frame_z.push_back(x[mol_indx[0][i-1]][2]);
                      }
              coor_data_vec1_x.push_back(per_frame_x);
              coor_data_vec1_y.push_back(per_frame_y);
              coor_data_vec1_z.push_back(per_frame_z);
              per_frame_x.clear();
              per_frame_y.clear();
              per_frame_z.clear();
              for(int i=0; i<mol_indx[1].size(); i++){
                      per_frame_x.push_back(x[mol_indx[1][i-1]][0]);
                      per_frame_y.push_back(x[mol_indx[1][i-1]][1]);
                      per_frame_z.push_back(x[mol_indx[1][i-1]][2]);
                      }
              coor_data_vec2_x.push_back(per_frame_x);
              coor_data_vec2_y.push_back(per_frame_y);
              coor_data_vec2_z.push_back(per_frame_z);
              time_last=time;
              nframe++;
              }
              }
     dt=(time_last-time_first)/(nframe-1);
     }
     
void calculate(int corn){
     omp_set_num_threads(corn);
     ofstream out_put(output_file.c_str(), ofstream::out);
     
     int ntime=floor((time_last-time_first)/dt);
     
     if(ntime<1000){
                    cerr<<"In order to ensure the accuracy of the results, enough trajectory should be given"<<endl;
                    if(ntime<1){
                                exit(0);
                                }
                    }
     
     for(int nt=0; nt<ntime; nt=nt+nspin){
             float sum_number=0.0;
             float sum1_number=0.0;
             float total_vec=0.0;
             float total_theta=0.0;
             for(int i=0; i<coor_data_vec1_x.size()-nt; i++){
                     #pragma omp parallel for reduction(+:total_theta) reduction(+:total_vec) reduction(+:sum1_number) reduction(+:sum_number)
                     for(int k=0; k<coor_data_vec1_x[i].size(); k++){
                             float vec1_x=coor_data_vec1_x[i][k]-coor_data_vec2_x[i][k];
                             float vec1_y=coor_data_vec1_y[i][k]-coor_data_vec2_y[i][k];
                             float vec1_z=coor_data_vec1_z[i][k]-coor_data_vec2_z[i][k];
                     
                             float vec2_x=coor_data_vec1_x[i+nt][k]-coor_data_vec2_x[i+nt][k];
                             float vec2_y=coor_data_vec1_y[i+nt][k]-coor_data_vec2_y[i+nt][k];
                             float vec2_z=coor_data_vec1_z[i+nt][k]-coor_data_vec2_z[i+nt][k];
                             
                             float value1=sqrt(vec1_x*vec1_x+vec1_y*vec1_y+vec1_z*vec1_z);
                             float value2=sqrt(vec2_x*vec2_x+vec2_y*vec2_y+vec2_z*vec2_z);
                             
                             if(value1!=0&&value2!=0){
                                                      total_theta +=(vec1_x*vec2_x+vec1_y*vec2_y+vec1_z*vec2_z)/value1/value2;
                                                      sum1_number++;
                                                      }
                             total_vec += vec1_x*vec2_x+vec1_y*vec2_y+vec1_z*vec2_z;
                             sum_number++;
                             }
                     }
             out_put<<(float)(nt+1)*dt/1000.0<<"\t"<<total_vec/sum_number<<"\t"<<total_theta/sum1_number<<endl;
             }
     
     }
     
int main(int argc, char *argv[]){
    int Max_core_number = sysconf(_SC_NPROCESSORS_CONF);
    
     if ((Threadnum <= 0) || (Threadnum > Max_core_number)){
                    Threadnum = Max_core_number;
                    } 
    Parse_Para(argc, argv);
    load_index(index_file.c_str());
    load_coor(traj_file);
    calculate(Threadnum);
    }
