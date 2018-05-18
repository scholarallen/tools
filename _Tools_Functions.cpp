/*
  Name: TOolFunctions
  Author: Allen Wang
  Date: 10/05/17 16:59
  Description: Functions:
                         BubbleSort --> 
                         DOT_Product -->
                         TMatrix -->
                         Matrix_MUL -->
                         Matrix_Sub_P -->
                         Matrix_Det -->
                         Matrix_Adjugate -->
                         Matrix_Inverse -->
*/

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
#include <algorithm>
#include <sys/stat.h>
#include <sys/types.h>

using namespace std;

/*******************************************************************************  
--- Sort array (array1) by bubble sort method, from smallest to largest. 
--- "rank1" is the original serial number.
--- "len" is the size of "array1" 
*******************************************************************************/

void BubbleSort(double *array1, int *rank1, int len){
     float temp;
     int temp_rank;
     
     for(int i = 0; i < len-1; i++){
             for(int j = 0; j < len-1-i; j++){
                     if(array1[j]>array1[j+1]){
                                               temp        = array1[j];
                                               array1[j]   = array1[j+1];
                                               array1[j+1] = temp;
                                               temp_rank   = rank1[j];
                                               rank1[j]    = rank1[j+1];
                                               rank1[j+1]  = temp_rank;
                                               }
                     }
             }
     }


/*******************************************************************************  
--- Dot production of two single-dimension deternminant. 
*******************************************************************************/

float DOT_Product(vector<float> Vector_A, vector<float> Vector_B, int Threadnum){
      if(Vector_A.size()!=Vector_B.size()){
                                           cout<<"the dimension of the two vector should be equal"<<endl;
                                           exit(0);
                                           }
      float Result=0.0;
      int Max_core_number = sysconf(_SC_NPROCESSORS_CONF);
      if ((Threadnum <= 0) || (Threadnum > Max_core_number)){
                     Threadnum = Max_core_number;
                     }
      omp_set_num_threads(Threadnum);
#pragma omp parallel for reduction(+:Result)
      for(int i=0; i<Vector_A.size(); i++){
              Result+=Vector_A[i]*Vector_B[i];
              }
      return Result;
      }


/*******************************************************************************  
--- Calculate the transposed matrix with two same dimensions. 
*******************************************************************************/
     
vector< vector< float > > TMatrix(vector< vector< float > > Matrix_A){
        vector< vector< float > > Matrix_B;
        for(int i=0; i<Matrix_A[0].size(); i++){
                vector<float> Per_Colum;
                for(int k=0; k<Matrix_A.size(); k++){
                        Per_Colum.push_back(Matrix_A[k][i]);
                        }
                Matrix_B.push_back(Per_Colum);
                }
        return Matrix_B;
        }

/*******************************************************************************  
--- Matrix multiplication. 
*******************************************************************************/

vector< vector< float > > Matrix_MUL(vector< vector< float > > Vector_A, vector< vector< float > > Vector_B, int Threadnum){
      
      vector< vector< float > > Vector_TB, Matrix_Total;
      Vector_TB=TMatrix(Vector_B);
      
      for(int i=0; i<Vector_A.size(); i++){
              vector<float> Per_Row;
              for(int k=0; k<Vector_TB.size(); k++){
                      Per_Row.push_back(DOT_Product(Vector_A[i],Vector_TB[k], Threadnum));
                      }
              Matrix_Total.push_back(Per_Row);
              }
      return Matrix_Total;
      }


/*******************************************************************************  
--- Calculation of Algebraic cofactor . 
*******************************************************************************/

vector< vector< float > > Matrix_Sub_P(vector< vector< float > > Vector_A, int nr, int nc){
        vector< vector< float > > Return_V;
        for(int i=0; i<Vector_A.size(); i++){
                if(i!=nr){
                          vector<float> Per_Line;
                          for(int k=0; k<Vector_A.size();k++){
                                  if(k!=nc){Per_Line.push_back(Vector_A[i][k]);}
                                  }
                          Return_V.push_back(Per_Line);
                          }
                }
        return Return_V;
        }
        
 
/*******************************************************************************  
--- Calculation of Algebraic cofactor, just used for "Matrix_Det". 
*******************************************************************************/
          
vector< vector< float > > Matrix_Sub(vector< vector< float > > Vector_A, int nr, int nc){
        vector< vector< float > > Return_V;
        for(int i=0; i<Vector_A.size(); i++){
                if(i!=nr){
                          vector<float> Per_Line;
                          for(int k=0; k<Vector_A.size();k++){
                                  int kk=pow(-1,nc);
                                  if(k!=nc){
                                            if(i==1){Per_Line.push_back(kk*Vector_A[0][nc]*Vector_A[i][k]);}else{
                                                                                                                 Per_Line.push_back(Vector_A[i][k]);
                                                                                                                 }
                                            }
                                  }
                          Return_V.push_back(Per_Line);
                          }
                }
        return Return_V;
        }


/*******************************************************************************  
--- Used the recursived method to calculate the determinant value. 
*******************************************************************************/
float Matrix_Det(vector< vector< float > > Vector_A){
      if(Vector_A.size()!=Vector_A[0].size()){
                                              cout<<"Errors::(Matrix_Det)The matrix must be a square"<<endl;
                                              exit(0);
                                              }
      if(Vector_A.size()==2){
                             return Vector_A[0][0]*Vector_A[1][1]-Vector_A[0][1]*Vector_A[1][0];
                             }else{
                                   float Result=0.0;
                                   vector< vector< float > > New_Matrix;
                                   for(int i=0; i<Vector_A.size(); i++){
                                           New_Matrix=Matrix_Sub(Vector_A, 0, i);
                                           Result+=Matrix_Det(New_Matrix);
                                           }
                                   return Result;
                                   }
      
      }
      

/*******************************************************************************  
--- Calculate the adjugate matrix. 
*******************************************************************************/

vector< vector< float > > Matrix_Adjugate(vector< vector< float > > Vector_A){
        if(Vector_A.size()!=Vector_A[0].size()){
                                                cout<<"Errors::(Matrix_Adjugate)The matrix must be a square"<<endl;
                                                exit(0);
                                                }
        vector< vector< float > > Return_Matrix;
        if(Vector_A.size()==2){
                               vector<float> Temp_Matrix1;
                               Temp_Matrix1.push_back(Vector_A[1][1]);
                               Temp_Matrix1.push_back(-1*Vector_A[0][1]);
                               Return_Matrix.push_back(Temp_Matrix1);
                               Temp_Matrix1.clear();
                               Temp_Matrix1.push_back(-1*Vector_A[1][0]);
                               Temp_Matrix1.push_back(Vector_A[0][0]);
                               Return_Matrix.push_back(Temp_Matrix1);
                               }else{
                                     for(int i=0; i<Vector_A.size(); i++){
                                             vector<float> Per_Line;
                                             for(int k=0; k<Vector_A[i].size(); k++){
                                                     int kk=pow(-1,i+k);
                                                     vector< vector< float > > Temp_Matrix;
                                                     Temp_Matrix=Matrix_Sub_P(Vector_A,k,i);
                                                     Per_Line.push_back(kk*Matrix_Det(Temp_Matrix));
                                                     }
                                             Return_Matrix.push_back(Per_Line);
                                             }
                                     }
        
        return Return_Matrix;
        }


/*******************************************************************************  
--- Calculate the inverse matrix. 
*******************************************************************************/

vector< vector< float > > Matrix_Inverse(vector< vector< float > > Vector_A){
        float Dot_Value=Matrix_Det(Vector_A);
        if(Dot_Value==0){
                         cerr<<"Error:: This matrix is not a inverse matrix."<<endl;
                         exit(0);
                         }
        vector< vector< float > > Vector_B;
        Vector_B=Matrix_Adjugate(Vector_A);
        for(int i=0; i<Vector_B.size(); i++){
                for(int k=0; k<Vector_B[i].size(); k++){
                        Vector_B[i][k]=Vector_B[i][k]/Dot_Value;
                        }
                }
        return Vector_B;
        }
