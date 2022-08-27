
#include "Self_Define_Functions.h"
#include<iostream>

using namespace std;

int main(int argc, char* argv[])
{
    int p, q, V, i, j, k,t;                                           //ѭ������
    int func_num;                                                //�������
    double final_val;                                           //������Ӧ��
    char fitness_file_name[256];                         //�ļ�����
    int run_index;                                               //���д���
    const int num_of_tiers = 7;                                   //������ÿ������Ⱥ�����ӵķֳɵĲ�����
    const int hierarchical = 4;                                     //�����������Ⱥ�ĸ�����
    int fitness_counter;                                      //��Ӧ�ȼ������
    int current_generation;                                //��ǰ����
    const double phi = 4.1;
    const double lambda = 1.2;
    const int subpopulation_size = Population_Size / hierarchical;//����Ⱥ����Ⱥ��С
    int subpopulation_length = subpopulation_size;            //����Ⱥ�е����Ӹ���������Ⱥ�ĳ��ȣ�
    const int exchange_frequency = 30;                                        //����Ƶ��
    const int reinitialization_frequency = 80;                                //���³�ʼ��Ƶ��
    double r1,r2;
    const double lowerbound=-100;                                                      //������Сֵ
    const double upperbound=100;                                                       //�������ֵ

    double** subpopulation_tier_gbest = new double* [num_of_tiers];//�洢��ǰ����Ⱥÿ��ÿ��subwarm��gbest
    double** population = new double* [Population_Size];//�洢������Ⱥ��ÿ�����ӵ�λ��������
    double** speed = new double* [Population_Size];//�洢ÿ�����ӵ��ٶ�����
    double** personal_best = new double* [Population_Size];//�洢ÿ�����ӵ�pbest��λ����Ϣ
    int** subpopulation = new int* [hierarchical];//�洢ÿ������Ⱥ��������population�е�λ�ã����Ϊ��i������Ⱥ�еĵ�j������Ϊt��tΪ��������population�е�λ�ã�
    double* Admission_threshold = new double[hierarchical];//ÿ������Ⱥ��׼����ֵ
    int** Buffer_fields = new int* [hierarchical];//ÿ������Ⱥ��׼�뻺����
    double* current_subpopulation_results = new double[subpopulation_size];//�洢��ǰ����Ⱥ��������Ӧ��
    double* temp_subpopulation_result = new double[subpopulation_size];//��ʱ�洢��ǰ����Ⱥ�����ӵ���Ӧ��                 
    double* temp_buffer_field_result = new double[Population_Size];//��ʱ�洢�����������ӵ���Ӧ��
    int* buffer_field_length = new int[Population_Size];//�洢�����������Ӹ���
    int*** subpopulation_gbest_index = new int** [hierarchical];//�洢ÿ������Ⱥ��ÿһ����ÿ��subwarms�����Ž�
    double* results = new double[Population_Size];//�洢�������ӵ���Ӧ��    
    double* personal_best_results = new double[Population_Size];//�洢ÿ�����ӵ���ʷ�����Ӧ��
    double* acceleration = new double[num_of_tiers];//�洢ÿ������Ⱥ�ļ���ϵ��
    double* temp1 = new double[Population_Size];
    double* temp2 = new double[Population_Size];
  
     //����ռ�
    for (i = 0; i < Population_Size; ++i)
    {
        population[i] = new double[dim];
        speed[i] = new double[dim];
        personal_best[i] = new double[dim];
    }

    for (i = 0; i < hierarchical; i++)
    {
        subpopulation[i] = new int[subpopulation_size];
        Buffer_fields[i] = new int[Population_Size];
    }

    for (i = 0; i < num_of_tiers; i++)
    {
        subpopulation_tier_gbest[i] = new double[subpopulation_size];
    }
    
    for (i = 0; i < hierarchical; i++)
    {
        subpopulation_gbest_index[i] = new int* [num_of_tiers];
        for (j = 0; j < num_of_tiers; j++)
        {
            subpopulation_gbest_index[i][j] = new int[subpopulation_size];
        }
    }

    //��ʼ��ÿ��ļ���ϵ��
    double temp = phi / (2 * lambda + num_of_tiers - 2);
    acceleration[0] = lambda * temp;                                                   //�����   
    acceleration[num_of_tiers - 1] = lambda * temp;                           //��num_of_tiers-1��
    for (i = 1; i < num_of_tiers - 1; i++)
    {
        acceleration[i] = temp;
    }
    boost::mt19937 generator(time(0) * rand());
    boost::uniform_real<> uniform_real_generate_r(0, 1);
    boost::variate_generator< boost::mt19937&, boost::uniform_real<> > random_real_num_r(generator, uniform_real_generate_r);//to generate a random number within [0,1]
   
    boost::uniform_real<> uniform_real_generate_x(lowerbound, upperbound);
    boost::variate_generator< boost::mt19937&, boost::uniform_real<> > random_real_num_x(generator, uniform_real_generate_x);

    //����ʵ��
    for (V = 0; V < 30; V++)
    {
       
        func_num = V + 1;
        sprintf_s(fitness_file_name, "F%d/convergency.txt", func_num);
        string FileName = string(fitness_file_name);
        ofstream  out_fitness(FileName);
        cout << "Start to optimize the function" << func_num << endl;
        for (run_index = 1; run_index <= 10; run_index++)
        {
           

            fitness_counter = 0;
            current_generation = 1;

            //��ʼ��population��speed
            for (i = 0; i < dim; ++i)
            {
                boost::uniform_real<> uniform_real_generate_x(lowerbound, upperbound);
                boost::variate_generator< boost::mt19937&, boost::uniform_real<> > random_real_num_x(generator, uniform_real_generate_x);

                for (j = 0; j < Population_Size; ++j)
                {
                    personal_best[j][i] = population[j][i] = random_real_num_x();
                    speed[j][i] = 0;
                }
            }


            //��������Ⱥ���з��飬�ֳ�H��subpopulation
            i = 0;
            j = 0;
            for (p = 0; p < Population_Size; p++)                   //p����������population�е�λ��
            {
                if (j == subpopulation_size - 1)                   //����ﵽsubpopulation_size������
                {
                    subpopulation[i][j] = p;
                    i++;
                    j = 0;
                }
                else
                {
                    subpopulation[i][j] = p;
                    j++;
                }
            }

            //��Ӧ�ȼ���
            for (i = 0; i < Population_Size; i++)
            {
                cec14_test_func(population[i], &results[i], dim, 1, func_num);
                results[i] = results[i] - (func_num) * 100;
                // cout<<results[i]<<endl;
            }

            //��ʼ��ÿ�����ӵ���ʷ������Ӧ��
            memcpy(personal_best_results, results, sizeof(double) * Population_Size);

            //Ѱ��������Ӧ��
            final_val = results[0];
            //cout << results[0]<<endl;

            for (i = 0; i < hierarchical; i++)
            {
                subpopulation_length = subpopulation_size;

                //��ʼ����0���������Ӧ���Լ�������Ӧ��������population�е�λ��

                for (j = 0; j < subpopulation_size; j++)
                {
                    subpopulation_tier_gbest[0][j] = results[subpopulation[i][j]];
                    subpopulation_gbest_index[i][0][j] = subpopulation[i][j];
                    //cout << results[subpopulation[i][j]]<<" ";
                }
                //�������{1��2��2��2��2��2��2}�ķֲ�������k���subwarm����������k-1���subwarm���ɣ���length��Ҫÿ�γ���2
                //��0���Ѿ�������ϣ��ʴ˴�Ҫ��2
                subpopulation_length /= 2;

                //Ѱ��ÿ��ÿ��subwarm��gbest
                for (k = 1; k < num_of_tiers; k++)                          //k�����k��
                {
                    j = 0;                                                                     //j�����j��subwarm
                    t = 0;                                                                     //t�����t��subwarm
                    while (j < subpopulation_length)
                    {
                        //��k���j��subwarm���ɵ�k-1���t��subwarm�͵�k-1���t+1��subwarm��ɣ�ȡ��Сֵ
                        subpopulation_tier_gbest[k][j] = min(subpopulation_tier_gbest[k - 1][t], subpopulation_tier_gbest[k - 1][t + 1]);
                        if (subpopulation_tier_gbest[k - 1][t] >= subpopulation_tier_gbest[k - 1][t + 1])
                        {
                            subpopulation_gbest_index[i][k][j] = subpopulation_gbest_index[i][k - 1][t + 1];
                        }
                        else  subpopulation_gbest_index[i][k][j] = subpopulation_gbest_index[i][k - 1][t];
                        //cout << gbest_index[i][k][j] << " ";
                        j++;
                        t += 2;
                    }
                    //cout << endl;
                    subpopulation_length /= 2;
                }

                //ÿһ��subpopulation�еĵ�k-1��ĵ�0��subwarm����һsubpopulation��ȫ������
                // cout << gbest[k - 1][0] << endl;
                if (final_val > subpopulation_tier_gbest[num_of_tiers - 1][0])
                {
                    final_val = subpopulation_tier_gbest[num_of_tiers - 1][0];
                    //gl_best = gbest_index[i][k - 1][0];
                }
                //cout << final_val << endl;
            }
            while (fitness_counter < MAX_FV)
            {


                if (current_generation % exchange_frequency == 0)
                {
                    for (i = 0; i < hierarchical; i++)
                    {
                        for (j = 0; j < subpopulation_size; j++)
                        {
                            current_subpopulation_results[j] = results[subpopulation[i][j]];
                        }
                        sort(current_subpopulation_results, current_subpopulation_results + subpopulation_size);
                        Admission_threshold[i] = current_subpopulation_results[subpopulation_size / 2];
                    }

                    for (i = hierarchical - 1; i >= 1; i--)
                    {
                        p = 0;
                        for (j = 0; j < i; j++)
                        {
                            for (t = 0; t < subpopulation_size; t++)
                            {
                                if (results[subpopulation[j][t]] < Admission_threshold[i])
                                {
                                    Buffer_fields[i][p] = subpopulation[j][t];
                                    //cout << Buffer_fields[i][p] << endl;
                                    p++;
                                }
                            }
                        }
                        buffer_field_length[i] = p;
                        // cout << sublength[i]<<endl;
                    }
                    for (i = hierarchical - 1; i >= 1; i--)
                    {

                        for (j = 0; j < buffer_field_length[i]; j++)
                        {
                            temp_buffer_field_result[j] = results[Buffer_fields[i][j]];
                            temp1[j] = results[Buffer_fields[i][j]];
                        }
                        for (j = 0; j < subpopulation_size; j++)
                        {
                            temp_subpopulation_result[j] = results[subpopulation[i][j]];
                            temp2[j] = results[subpopulation[i][j]];
                        }
                        replaceinfo(results, Buffer_fields, subpopulation[i], i, temp_subpopulation_result, temp_buffer_field_result, subpopulation_size, buffer_field_length[i], temp1,temp2);

                    }

                    if (current_generation % reinitialization_frequency == 0)
                    {
                        for (i = 0; i < dim; ++i)
                        {
                            for (j = 0; j < subpopulation_size; ++j)
                            {
                                population[subpopulation[0][j]][i] = random_real_num_x();
                            }
                        }
                        for (j = 0; j < subpopulation_size; j++)
                        {
                            cec14_test_func(population[subpopulation[0][j]], &results[subpopulation[0][j]], dim, 1, func_num);
                            results[subpopulation[0][j]] = results[subpopulation[0][j]] - (func_num) * 100;
                            if (results[subpopulation[0][j]] < personal_best_results[subpopulation[0][j]])
                            {
                                personal_best_results[subpopulation[0][j]] = results[subpopulation[0][j]];
                                memcpy(personal_best[subpopulation[0][j]], population[subpopulation[0][j]], sizeof(double) * dim);

                            }
                        }
                        fitness_counter += subpopulation_size;
                    }

                }


                //���¹�ʽ
                for (i = 0; i < hierarchical; i++)
                {
                    
                    for (j = 0; j < subpopulation_size; j++)
                    {
                        for (k = 0; k < dim; k++)
                        {
                            r1 = random_real_num_r();
                            r2 = random_real_num_r();
                            speed[subpopulation[i][j]][k] = weight * r1 * speed[subpopulation[i][j]][k] + acceleration[0] * r2 * (personal_best[subpopulation[i][j]][k] - population[subpopulation[i][j]][k]) + cal_glbest(j, subpopulation_tier_gbest, acceleration, num_of_tiers, population, subpopulation, i, k, subpopulation_gbest_index);
                        
                            //speed[subpopulation[i][j]][k] = weight * r * speed[subpopulation[i][j]][k] + c[0] * r * (pbest[subpopulation[i][j]][k] - population[subpopulation[i][j]][k]) + c[i] * r * (population[gl_best][k] - population[subpopulation[i][j]][k]);
                            population[subpopulation[i][j]][k] += speed[subpopulation[i][j]][k];
                            
                            if (population[subpopulation[i][j]][k] > upperbound)
                                population[subpopulation[i][j]][k] = upperbound;
                            else if (population[subpopulation[i][j]][k] < lowerbound)
                                population[subpopulation[i][j]][k] = lowerbound;
                               
                        }
                    }
                }


                //��Ӧ�ȼ���
                for (i = 0; i < Population_Size; i++)
                {
                    cec14_test_func(population[i], &results[i], dim, 1, func_num);
                    results[i] = results[i] - (func_num) * 100;
                    if (results[i] < personal_best_results[i])
                    {
                        personal_best_results[i] = results[i];
                        memcpy(personal_best[i], population[i], sizeof(double) * dim);

                    }
                }


                //����final_val
                for (i = 0; i < hierarchical; i++)
                {
                    subpopulation_length = subpopulation_size;

                    //��ʼ����0���������Ӧ���Լ�������Ӧ��������population�е�λ��

                    for (j = 0; j < subpopulation_size; j++)
                    {
                        subpopulation_tier_gbest[0][j] = results[subpopulation[i][j]];
                        subpopulation_gbest_index[i][0][j] = subpopulation[i][j];
                        //cout << results[subpopulation[i][j]]<<" ";
                    }
                    //�������{1��2��2��2��2��2��2}�ķֲ�������k���subwarm����������k-1���subwarm���ɣ���length��Ҫÿ�γ���2
                    //��0���Ѿ�������ϣ��ʴ˴�Ҫ��2
                    subpopulation_length /= 2;

                    //Ѱ��ÿ��ÿ��subwarm��gbest
                    for (k = 1; k < num_of_tiers; k++)                          //k�����k��
                    {
                        j = 0;                                                                     //j�����j��subwarm
                        t = 0;                                                                     //t�����t��subwarm
                        while (j < subpopulation_length)
                        {
                            //��k���j��subwarm���ɵ�k-1���t��subwarm�͵�k-1���t+1��subwarm��ɣ�ȡ��Сֵ
                            subpopulation_tier_gbest[k][j] = min(subpopulation_tier_gbest[k - 1][t], subpopulation_tier_gbest[k - 1][t + 1]);
                            if (subpopulation_tier_gbest[k - 1][t] >= subpopulation_tier_gbest[k - 1][t + 1])
                            {
                                subpopulation_gbest_index[i][k][j] = subpopulation_gbest_index[i][k - 1][t + 1];
                            }
                            else  subpopulation_gbest_index[i][k][j] = subpopulation_gbest_index[i][k - 1][t];
                            //cout << gbest_index[i][k][j] << " ";
                            j++;
                            t += 2;
                        }
                        //cout << endl;
                        subpopulation_length /= 2;
                    }
                    //ÿһ��subpopulation�еĵ�k-1��ĵ�0��subwarm����һsubpopulation��ȫ������
                    // cout << gbest[k - 1][0] << endl;
                    if (final_val > subpopulation_tier_gbest[k - 1][0])
                    {
                        final_val = subpopulation_tier_gbest[k - 1][0];
                        //gl_best = gbest_index[i][k - 1][0];
                    }
                    //cout << final_val << endl;
                }
                
                current_generation += 1;
                fitness_counter += Population_Size;
                
            }
            //out_fitness << final_val << endl;
            cout << final_val << endl;
        }
       
        
        cout << "The optimization of the function "<<func_num<<"  is finished!!" << endl;
       
        
    }
    
    for (i = 0; i < Population_Size; ++i)
    {
        delete[]population[i];
        delete[]speed[i];
        delete[]personal_best[i];

    }
    //cout << 1 <<endl;
    for (i = 0; i < hierarchical; i++)
    {
        delete[]Buffer_fields[i];
    }
    //cout << 2 << endl;
    delete[]population;
    delete[]speed;
    delete[]personal_best;
    //cout << 3 << endl;
    delete[]buffer_field_length;
    delete[]Buffer_fields;
    delete[]personal_best_results;
    //cout << 4 << endl;
    delete[]temp_subpopulation_result;
    delete[]temp_buffer_field_result;
    //cout << 5 << endl;
    delete[]current_subpopulation_results;
    //cout << 5.1 << endl;
    delete[]temp1;
    //cout << 6 << endl;
    for (i = 0; i < hierarchical; i++)
    {
        for (j = 0; j < num_of_tiers; j++)
        {
            delete[]subpopulation_gbest_index[i][j];
        }
        delete[]subpopulation_gbest_index[i];
    }
   // cout << 7 << endl;
    delete[]subpopulation_gbest_index;
    //cout << 8 << endl;
    for (i = 0; i < num_of_tiers; i++)
    {
        delete[]subpopulation_tier_gbest[i];
    }
    delete[]subpopulation_tier_gbest;
    //cout << 9 << endl;



    return 0;
}//end of main



