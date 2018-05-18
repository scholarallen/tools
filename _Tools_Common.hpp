/*
 * Name: Additional Tools
 * Copyright: 
 * Author: Allen Wang
 * Date: 28/12/16 09:38
 * Description: 
 *              Get_Random: generate random numbers (0~1). 
 *                          float: single number. 
 *                          vector: users can formulate the numbers. "numbs": The numbers of random numbers that you need.
 */
#ifndef _Common_Tools_H
#define _Common_Tools_H

#include <iostream>
#include <stdlib.h>
#include <fstream>
#include <sstream>
#include <stdio.h>
#include <string>
#include <vector>
#include <time.h>
#include <map>
#include <set>

#define random(x) (rand()%x)

using namespace std;

class _Tools_Additional{
      public:
      float Get_Random(){srand((int)time(0)); return random(1000)/1000.0;};
      vector <float> Get_Random(int numbs);
}; 
//

vector <float> _Tools_Additional::Get_Random(int numbs){
       srand((int)time(0));
       vector <float> Value_Random; 
       for(int i=0; i<numbs; i++)
               Value_Random.push_back(random(1000)/1000.0);
               
       return Value_Random;
       }


#endif /* _Common_Tools_H */
