
#include "Self_Define_Functions.h"
#include <math.h>
#include <iostream>
#include <fstream>
#include<algorithm>
using namespace std;



double cal_glbest(int index, double** gbest, double* c, int N, double** population,int **subpopulation,int h,int t,int ***gbest_index)
{
    boost::mt19937 generator(time(0) * rand());
    boost::uniform_real<> uniform_real_generate_r(0, 1);
    boost::variate_generator< boost::mt19937&, boost::uniform_real<> > random_real_num_r(generator, uniform_real_generate_r);
    
    double value = 0.0;
    int i, j, k,x;
    double r1;
    //r1 = random_real_num_r();
    x = subpopulation[h][index];
    for (k = 1; k < N; k++)
    {
        r1 = random_real_num_r();
        //if (index == 0 )  index = 0;
        if (index % 2 == 1)
        {
            index = (index - 1) / 2;
        }
        else index /= 2;
        //printf("%E\n", population[gbest_index[H][k][index]][t]);
        value += c[k] * r1 * (population[gbest_index[h][k][index]][t] - population[x][t]);
       
    }
    
    //value += c[1] * r * (population[gl_best][t] - population[subpopulation[h][index]][t]);
    return value;
}

void replaceinfo(double* results, int** Buffer, int *subpopulation, int hierarchical,double *subpopulation_fitness,double *buffer_field_fitness,int subpopulationsize,int length,double *temp_buffer_field_result,double *fitness2)
{
    int i, j,index1,index2;
    double temp;
    

    //将缓冲区的粒子进行排序
    //sort(buffer_field_fitness, buffer_field_fitness + length);
    
    //找到亚种群中适应度最大的一个粒子，将适应值赋值给temp
    
    temp = search_max_num(subpopulation_fitness, subpopulationsize);
    i = 0;
    while (buffer_field_fitness[i] < temp && i<length)
    {
        index1 = search_index(temp, subpopulation_fitness, subpopulationsize);  //寻找temp在亚种群中的位置
        index2 = search_index(buffer_field_fitness[i], temp_buffer_field_result, length);  //寻找当前缓冲区粒子在缓冲区的位置
        if (index1 == -1)
        {
            break;
        }
        else
        {
            subpopulation[index1] = Buffer[hierarchical][index2];
            subpopulation_fitness[index1] = Buffer[hierarchical][index2];
            i++;
            temp = search_max_num(subpopulation_fitness, subpopulationsize);
        }
    }
    
    
    
    
    
    /*
    sort(buffer_field_fitness, buffer_field_fitness + length);
    sort(subpopulation_fitness, subpopulation_fitness + subpopulationsize);
    for (i = 0,j=subpopulationsize-1; i < subpopulationsize / 2 && i<length; i++,j--)
    {
        index1 = search_index(subpopulation_fitness[j], fitness2, subpopulationsize);
        index2 = search_index(buffer_field_fitness[i], fitness1, length);
        subpopulation[index1] = Buffer[hierarchical][index2];
    }
    */
    
}


int search_index(double x, double* fitness,int length)
{
    int i;
    for (i = 0; i < length; i++)
    {
        if (x == fitness[i])       return i;
    }
    return -1;
}

int search_max_num(double *goal,int length)
{
    sort(goal, goal + length);
    return goal[length - 1];
}






