#include <stdio.h>
#include "GA.h"
#include "TridiagonalSystemsOfEquations.h"
#include <time.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string>

using namespace std;

double x_matrix[n + 1];
//�Գ����ԽǷ��̵��������b
double b_matrix[n - 1];
//�Գ����ԽǷ��̵��������A
double A_matrix[3 * n - 5];
double u_observe[n - 1];


int main()
{
    //int random_num = (unsigned)time(NULL);
    int random_num = 1;
    srand(random_num);
    //��ʼ��4������
    for (int i = 0; i <= n; i++)
    {
        x_matrix[i] = i * x_gap;
    }
    for (int i = 0; i <= n - 2; i++)
    {
        u_observe[i] = (1 + delta * ((rand() % noise_point_num) * noise_gap - 1)) * u(x_matrix[i + 1]);
        //cout << u_observe[i]/u(x_matrix[i + 1]) << " ";
    }
    for (int i = 0; i <= n - 2; i++)
    {
        b_matrix[i] = 2 * x_gap * x_gap * f(x_matrix[i + 1]);
    }
    for (int i = 0; i < 3 * n - 5; i++)
    {
        A_matrix[i] = 0;
    }
    int LoopNum =  300000;
    int LoopNum2 = 300000;
    bool BreakLoop;
    clock_t start, end;
    int loop_num_real;
    int test_num = 10;
    string str = "result";

    ofstream fout_statis;
    fout_statis.open("statistics.txt");
    fout_statis << "����" << " " << "1.1REQ" << " " << "1.1REU" << " " << "1.1time" << " " << "1.2REQ" << " " << "1.2REU" << " " << "1.2time" << " " << "1.3REQ" << " " << "1.3REU" << " " << "1.3time" << " " << "2.1REQ" << " " << "2.1REU" << " " << "2.1time" << " " << "2.2REQ" << " " << "2.2REU" << " " << "2.2time" <<" " << "2.3REQ" << " " << "2.3REU" << " " << "2.3time" << endl;

    for (int test_num_real = 0; test_num_real < test_num; test_num_real++)
    {
        cout << "this is the " << test_num_real << " loop" << endl;

        cout << "��һʵ�鿪ʼ" << endl;
        srand(random_num + test_num_real);
        GenGroup* Group = new GenGroup;
        loop_num_real = 0;
        start = clock();
        //��Ҫ��������ʱ��ĳ����
        BreakLoop = Group->CheckStop();
        while (loop_num_real < 10000 && !BreakLoop)
        {
            loop_num_real += 1;
            BreakLoop = Group->evolution();
            //BreakLoop = Group->variation();
        }
        end = clock();
        /*
        if (Group->REU<0.01)//��Ϊ��������ֲ�����
        {
            continue;
        }
        
        */

        fout_statis << test_num_real << " ";
        cout << "��һ�غ����չʾ��" << endl;
        cout << "���лغ�����" << loop_num_real << endl;
        cout << "REQ:" << Group->REQ << "REU:" << Group->REU << endl;
        cout << "��һ�غ�����ʱ��" << (double)(end - start) / CLOCKS_PER_SEC << endl;
        fout_statis << Group->REQ << " " << Group->REU << " ";
        fout_statis << (double)(end - start) / CLOCKS_PER_SEC << " ";


        while (loop_num_real < 200000 && !BreakLoop)
        {
            loop_num_real += 1;
            BreakLoop = Group->evolution();
            //BreakLoop = Group->variation();
        }
        end = clock();
        cout << "�ڶ��غ����չʾ��" << endl;
        cout << "���лغ�����" << loop_num_real << endl;
        cout << "REQ:" << Group->REQ << "REU:" << Group->REU << endl;
        cout << "��һ�ӵڶ��غ�����ʱ��" << (double)(end - start) / CLOCKS_PER_SEC << endl;
        fout_statis << Group->REQ << " " << Group->REU << " ";
        fout_statis << (double)(end - start) / CLOCKS_PER_SEC << " ";


        while (loop_num_real < LoopNum && !BreakLoop)
        {
            loop_num_real += 1;
            BreakLoop = Group->evolution();
            //BreakLoop = Group->variation();
        }
        end = clock();
        cout << "�����غ����չʾ��" << endl;
        cout << "���лغ�����" << loop_num_real << endl;
        cout << "REQ:" << Group->REQ << "REU:" << Group->REU << endl;
        cout << "��һ�ӵڶ��ӵ����غ�����ʱ��" << (double)(end - start) / CLOCKS_PER_SEC << endl;
        fout_statis << Group->REQ << " " << Group->REU << " ";
        fout_statis << (double)(end - start) / CLOCKS_PER_SEC << " ";

        cout << "output our result and the explicit q value" << endl;
        cout << "��10���ľ�ϸ��" << endl;
        const int x_output_len = 10 * n + 1;
        double x_output_gap = x_gap / 10;
        double x_matrix_output[x_output_len];
        for (int i = 0; i < x_output_len; i++)
        {
            x_matrix_output[i] = i * x_output_gap;
        }

        //�����������txt��
        ofstream fout;


        fout.open(str + to_string(test_num_real) + ".1.txt");

        for (int i = 0; i < x_output_len; i++)
        {
            fout << x_matrix_output[i] << " ";
        }
        fout << endl;
        for (int i = 0; i < x_output_len; i++)
        {
            fout << Group->Gens[0].q_value(x_matrix_output[i]) << " ";
        }
        fout << endl;
        for (int i = 0; i < x_output_len; i++)
        {
            fout << q(x_matrix_output[i]) << " ";
        }
        fout << endl;
        fout.close();//��ɺ󣬹ر�TXT�ļ�
        //��һ��ʵ�����
        delete Group;




        //�ڶ���ʵ�飬ʹ����ͬ����
        cout << "�ڶ�ʵ�鿪ʼ" << endl;
        srand(random_num + test_num_real);
        Group = new GenGroup;
        BreakLoop = Group->CheckStop();
        start = clock();
        loop_num_real = 0;
        while (loop_num_real < 10000 && !BreakLoop)
        {
            loop_num_real += 1;
            BreakLoop = Group->evolution();
            //BreakLoop = Group->variation();
        }
        end = clock();
        cout << "��һ�غ����չʾ��" << endl;
        cout << "���лغ�����" << loop_num_real << endl;
        cout << "REQ:" << Group->REQ << "REU:" << Group->REU << endl;
        cout << "��һ�غ�����ʱ��" << (double)(end - start) / CLOCKS_PER_SEC << endl;
        fout_statis << Group->REQ << " " << Group->REU << " ";
        fout_statis << (double)(end - start) / CLOCKS_PER_SEC << " ";


        while (loop_num_real < 200000 && !BreakLoop)
        {
            loop_num_real += 1;
            BreakLoop = Group->evolution();
            BreakLoop = Group->variation();
        }
        end = clock();
        cout << "�ڶ��غ����չʾ��" << endl;
        cout << "���лغ�����" << loop_num_real << endl;
        cout << "REQ:" << Group->REQ << "REU:" << Group->REU << endl;
        cout << "��һ�ӵڶ��غ�����ʱ��" << (double)(end - start) / CLOCKS_PER_SEC << endl;
        fout_statis << Group->REQ << " " << Group->REU << " ";
        fout_statis << (double)(end - start) / CLOCKS_PER_SEC << " ";

        while (loop_num_real < LoopNum2 && !BreakLoop)
        {
            loop_num_real += 1;
            BreakLoop = Group->evolution();
            BreakLoop = Group->variation();
        }
        end = clock();
        cout << "�����غ����չʾ��" << endl;
        cout << "���лغ�����" << loop_num_real << endl;
        cout << "REQ:" << Group->REQ << "REU:" << Group->REU << endl;
        cout << "��һ�ӵڶ��ӵ����غ�����ʱ��" << (double)(end - start) / CLOCKS_PER_SEC << endl;
        fout_statis << Group->REQ << " " << Group->REU << " ";
        fout_statis << (double)(end - start) / CLOCKS_PER_SEC << endl;

        cout << "output our result and the explicit q value" << endl;
        cout << "��10���ľ�ϸ��" << endl;
        for (int i = 0; i < x_output_len; i++)
        {
            x_matrix_output[i] = i * x_output_gap;
        }
        //�����������txt��
        fout.open(str + to_string(test_num_real) + ".2.txt");

        for (int i = 0; i < x_output_len; i++)
        {
            fout << x_matrix_output[i] << " ";
        }
        fout << endl;
        for (int i = 0; i < x_output_len; i++)
        {
            fout << Group->Gens[0].q_value(x_matrix_output[i]) << " ";
        }
        fout << endl;
        for (int i = 0; i < x_output_len; i++)
        {
            fout << q(x_matrix_output[i]) << " ";
        }
        fout << endl;
        fout.close();//��ɺ󣬹ر�TXT�ļ�
        delete Group;

        cout << endl;
    }
    fout_statis <<"delta" <<" "<< "beta" << " " <<"EliteNumber"<<" "<<"ChildrenNumber" << " " << "variation_rate" << " " <<"variation_gap_num_scale_rate" << endl;
    fout_statis << delta << " " << beta << " " << EliteNumber  << " "<< ChildrenNumber << " " <<  variation_rate  << " " << variation_gap_num_scale_rate<< endl;
    return 0;
}