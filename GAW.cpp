/*
 * The gentic algorithm operator.
 * by Honglei Wang, 2019.10.11
 * E-Mail :: scholar_allen@yeah.net
 */
 
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <map>
#include <cmath>
#include <string>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>

#define random(x) (rand()%x)
#define Dlength 100000
#define ND_ratio 0.7
#define ND_Error 0.6

using namespace std;

float DDratio=0.1;
string intTstring(int input_int){
	ostringstream os;
	os<<input_int;
	return os.str();
}

int Exchang_Genemone(float ratio, float& Fvec, float& Mvec, float minV, float maxV){
	
	if(ratio>1){
		cerr<<"ERROR:: The ratio of mutation must be less than or equal to 1.0"<<endl;
		return 0;
	}

	int Fvec1=floor(Fvec);
	int Fvec2=floor((Fvec-Fvec1)*Dlength);
	int Mvec1=floor(Mvec);
	int Mvec2=floor((Mvec-Mvec1)*Dlength);
	
	string sFvec1=intTstring(Fvec1);
	string sFvec2=intTstring(Fvec2);
	string sMvec1=intTstring(Mvec1);
	string sMvec2=intTstring(Mvec2);
	
	int Real_Dimension=(sFvec1.length()>sMvec1.length())?sFvec1.length():sMvec1.length();
	
	if(Real_Dimension-sFvec1.length()>0){
		string Temp_Char;
		for(int i=0; i<Real_Dimension-sFvec1.length(); i++){
			Temp_Char+="0";
		}
		sFvec1=Temp_Char+sFvec1;
	}
	
	if(Real_Dimension-sMvec1.length()>0){
		string Temp_Char;
		for(int i=0; i<Real_Dimension-sMvec1.length(); i++){
			Temp_Char+="0";
		}
		sMvec1=Temp_Char+sMvec1;
	}
	
	for(int i=0; i<Real_Dimension; i++){
		float real_ratio=random(1000)/1000.0;
		if(real_ratio<ratio){
			char Temp_Gene=sMvec1[i];
			sMvec1[i]=sFvec1[i];
			sFvec1[i]=Temp_Gene;
		}
	}
	
	for(int i=0; i<sFvec2.length(); i++){
		float real_ratio=random(1000)/1000.0;
		if(real_ratio<ratio){
			char Temp_Gene=sMvec2[i];
			sMvec2[i]=sFvec2[i];
			sFvec2[i]=Temp_Gene;
		}
	}
	
	Fvec=atof(sFvec1.c_str())+atof(sFvec2.c_str())/Dlength;
	Mvec=atof(sMvec1.c_str())+atof(sMvec2.c_str())/Dlength;
	if(Fvec>maxV){Fvec=minV+(maxV-minV)*random(1000)/1000.0;}
	if(Fvec<minV){Fvec=minV+(maxV-minV)*random(1000)/1000.0;}
	if(Mvec>maxV){Mvec=minV+(maxV-minV)*random(1000)/1000.0;}
	if(Mvec<minV){Mvec=minV+(maxV-minV)*random(1000)/1000.0;}
	return 0;
}

int Mutation_Genemone(float ratio, float& Fvec, float minV, float maxV){
/*the value of mutation defer to the normal distribution*/	
	if(ratio>1){
		cerr<<"ERROR:: The ratio of mutation must be less than or equal to 1.0"<<endl;
		return 0;
	}
	int ts=1;
	float real_ratio=random(1000)/1000.0;
	float ts_ratio=random(1000)/1000.0;
	if(ts_ratio>0.5){
		ts=-1;
	}
	
	if(minV>0){ts=1;}
	if(real_ratio<ratio){
		float Real_Error=log(ND_ratio+(1-ND_ratio)*random(1000)/1000.0);
		Fvec=ts*sqrt(-2*ND_Error*ND_Error*Real_Error)+Fvec;
	}
	if(Fvec>maxV){Fvec=minV+(maxV-minV)*random(1000)/1000.0;}
	if(Fvec<minV){Fvec=minV+(maxV-minV)*random(1000)/1000.0;}
	
	return 0;
}

int Initial(){
	srand((int)time(0)); 
}

/*******************************************************************************  
--- Sort array (array1) by bubble sort method, from largest to smallest. 
*******************************************************************************/

int BubbleSort(vector<float>& Fitness, vector<float>& POPVec){
     
    int Index_lastExchange = 0;
    int sortBorder = Fitness.size();
	
	if(Fitness.size()!=POPVec.size()){
		cerr<<"ERROR: The length of vector(Fitness) must be equal to the length of POP "<<endl;
		return 0;
	}
	
    for(int i=0; i<sortBorder+1; i++){
        bool isSorted = true;
        for(int k=1; k<sortBorder; k++){
            float temp, rank_temp;
            if(Fitness[k]>Fitness[k-1]){
                temp = Fitness[k];
				rank_temp=POPVec[k];
                Fitness[k] = Fitness[k-1];
                Fitness[k-1] = temp;
				POPVec[k]=POPVec[k-1];
				POPVec[k-1]=rank_temp;
                isSorted = false;
                Index_lastExchange = k;
            }
        }
        sortBorder = Index_lastExchange;
        if(isSorted) {break;}
    }
	return 0;
}

int GAW(vector<float>& RFVec, vector<float> FITVec, float Ratio_EX, float Ratio_MT,float minV, float maxV){
/***RFVec : Parent
 ***FITVec: Fitness
 */
	BubbleSort(FITVec, RFVec);
	
	float sum_fitness=0.0;
	for(int i=0; i<FITVec.size(); i++){
		sum_fitness+=FITVec[i];
	}
	
	vector<float> Next_RFvec;
	for(int i=0; i<RFVec.size(); i++){
		Next_RFvec.push_back(RFVec[i]);
	}
	
	
	for(int isnPOP=1;isnPOP<RFVec.size();isnPOP++){
		int s1;
		vector<float> Next_RFvec1;
		for(int i=0; i<RFVec.size(); i++){
			Next_RFvec1.push_back(RFVec[i]);
		}
		float Total_Fitness=0;
		float s1_ratio=random(1000)/1000.0;
		int is=0;
		while(true){
			int th=random(FITVec.size()-1);
			if(th!=isnPOP){
				Total_Fitness+=FITVec[th]/sum_fitness;
			}
			if(Total_Fitness>s1_ratio){
				s1=th;
				break;
			}
			if(is>FITVec.size()){
				if(th!=isnPOP&&th>-1){
					s1=th;
					break;
				}
			}
			is++;
		}
		Exchang_Genemone(Ratio_EX, Next_RFvec1[s1], Next_RFvec1[isnPOP], minV, maxV);
		Mutation_Genemone(Ratio_MT,Next_RFvec1[isnPOP], minV, maxV);
		Next_RFvec[isnPOP]=Next_RFvec1[isnPOP];
	}

	RFVec=Next_RFvec;
	
}
