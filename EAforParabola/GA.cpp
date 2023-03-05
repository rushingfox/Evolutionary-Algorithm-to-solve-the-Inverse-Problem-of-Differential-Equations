#include "GA.h"
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <algorithm>
#include "TridiagonalSystemsOfEquations.h"
#include <iostream>
#include <vector>
using namespace std;

double beta = pow(10, -4);
double variation_rate = 0.005;
double variation_gap_num_scale_rate = 0.01;
int variation_gap_num_scale = alpha_max_num * variation_gap_num_scale_rate;
int variation_num = variation_rate * GroupScale;

//欧几里得距离计算函数
//欧式主函数
double EuclidDistance(double * m1,double * m2, int size)
{
    double result = 0;
    for (int i = 0; i < size; i++)
    {
        result += pow((m1[i] - m2[i]), 2);
    }
    result = pow(result, 0.5);
    return result;
}
double u(double x, double t)
{
    return sin(M_PI * x) * exp(t);
}
double q(double x)
{
    return exp(x);
}
//欧式函数重载——可以计算单向量和零向量的距离
double EuclidDistance(double* m1,int size)
{
    double result = 0;
    for (int i = 0; i < size; i++)
    {
        result += pow(m1[i], 2);
    }
    result = pow(result, 0.5);
    return result;
}


//刷新A矩阵
void A_matrix_change(Gen* g)
{
    for (int i = 0; i < A_matrix_len; i++)
    {
        int j = i / 3;
        if (i % 3 == 0)
        {
            A_matrix[i] = g->q[j] + 2 * g->q[j + 1] + g->q[j + 2];
        }
        else
        {
            A_matrix[i] = -1 * (g->q[j + 1] + g->q[j + 2]);
        }
    }
}


//根据求导规则求出来的右端函数
double f(double x, double t)
{
    return exp(t) * (sin(M_PI * x) - M_PI * exp(x) * cos(M_PI * x) + M_PI * M_PI * exp(x) * sin(M_PI * x));
}

//神奇，这个改进居然不work
/*

double fitness(Gen* g)
{
    double J = 0;
    //double* d;
    A_matrix_change(g);
    //double A_matrix_copy[A_matrix_len];
    //memcpy(A_matrix_copy, A_matrix, A_matrix_len * sizeof(double));
    double b_matrix_copy[n - 1];
    memcpy(b_matrix_copy, b_matrix, b_and_u_matrix_len * sizeof(double));
    //d = trde(b_and_u_matrix_len, A_matrix, b_matrix_copy);
    trde(b_and_u_matrix_len, A_matrix, b_matrix_copy);
    memcpy(g->u, b_matrix_copy, b_and_u_matrix_len * sizeof(double));
    J = pow(EuclidDistance(u_observe, g->u, b_and_u_matrix_len),2);
    J = J * x_gap;
    double P = 0;
    for (int i = 0; i < x_matrix_len-1; i++)
    {
        P += pow((g->q[i] - g->q[i + 1]), 2) / x_gap;
    }
    //cout << J + beta * P << "=J " << J << "+beta*P " << beta * P << endl;
    return J + beta * P;
}


*/


double fitness(Gen* g)
{
    double J = 0;
    double* d;
    A_matrix_change(g);
    double A_matrix_copy[A_matrix_len];
    memcpy(A_matrix_copy, A_matrix, A_matrix_len * sizeof(double));
    double b_matrix_copy[n - 1];
    memcpy(b_matrix_copy, b_matrix, b_and_u_matrix_len * sizeof(double));
    d = trde(b_and_u_matrix_len, A_matrix_copy, b_matrix_copy);
    memcpy(g->u, d, b_and_u_matrix_len * sizeof(double));
    J = pow(EuclidDistance(u_observe, g->u, b_and_u_matrix_len), 2);
    J = J * x_gap;
    double P = 0;
    for (int i = 0; i < x_matrix_len - 1; i++)
    {
        P += pow((g->q[i] - g->q[i + 1]), 2) / x_gap;
    }
    //cout << J + beta * P << "=J " << J << "+beta*P " << beta * P << endl;
    return J + beta * P;
}


double Gen::q_value(double x)
{
    int i = (long)((x - x_min) / x_gap);
    return q[i] * (x_min + (i + 1) * x_gap - x) / x_gap + q[i + 1] * (x - x_min - i * x_gap) / x_gap;
}

double Gen::GenUpdate()
{
    GenFitness = fitness(this);
    return GenFitness;
}

void Gen::GenRandomGenerate()
{
    for (int i = 0; i < x_matrix_len; i++)
    {
        this->q[i] = (rand() % (alpha_max_num + 1)) * alpha_gap + alpha_low;
    }
    this->GenUpdate();
}

bool Gen::operator==(const Gen& g)
{
    for (int i = 0; i <x_matrix_len ; i++)
    {
        if (this->q[i]!=g.q[i])
        {
            return false;
        }
    }
    return true;
}

void Gen::SmoothQ()
{
    for (int i = 1; i <= n-1; i++)
    {
        this->q[i] = (this->q[i - 1] + this->q[i] + this->q[i + 1]) / 3;
    }
}

void Gen::GenPrint()
{
    cout << "the q is:" << endl;
    for (int i = 0; i < x_matrix_len; i++)
    {
        cout << this->q[i] << " ";
    }
    cout << endl;
}


bool Gen::operator< (const Gen & s)
{
    return this->GenFitness < s.GenFitness;
}


//其实就是先拷贝q再生成其他的
void Gen::GenCopy(const Gen& g_origin)
{
    for (int i = 0; i < x_matrix_len; i++)
    {
        this->q[i] = g_origin.q[i];
    }
    for (int i = 0; i < b_and_u_matrix_len; i++)
    {
        this->u[i] = g_origin.u[i];
    }
    this->GenFitness = g_origin.GenFitness;
}

GenGroup::GenGroup()
{
    for (int i = 0; i < GroupScale; i++)
    {
        Gens[i].GenRandomGenerate();
    }
    GensUpdate();
}

void GenGroup::GensUpdate()
{
    /*
    for (int i = 0; i < GroupScale; i++)
    {
        this->Gens[i].GenUpdate();
    }
    */

    sort(this->Gens, this->Gens + GroupScale);
    //cout << Gens[0].GenFitness << " " << Gens[GroupScale - 1].GenFitness;
    REQ = EuclidDistance(this->Gens[0].q, this->Gens[GroupScale - 1].q, x_matrix_len) / EuclidDistance(this->Gens[0].q, x_matrix_len);
    REU = EuclidDistance(this->Gens[0].u, u_observe, b_and_u_matrix_len) / EuclidDistance(u_observe, b_and_u_matrix_len);
    //cout << " REQ:" << REQ << " REU:" << REU << endl;
}

void GenGroup::ChildrenSpaceGenerate(double * p)
{
    double sum = 0;
    for (int i = 0; i < dimension-1; i++)
    {
        double temp ;
        do
        {
            temp = (rand() % a_n) * a_gap + a_low;
        } while (!((temp + sum)>=(1-a_up)&&(temp + sum)<=(1-a_low)));
        p[i] = temp;
        sum += p[i];
        //cout << p[i] << " ";
    }
    p[dimension - 1] = 1 - sum;
    //cout << p[dimension - 1];
}

bool GenGroup::CheckStop()
{
    if (REQ <= varepsilonForREQ && REU <= varepsilonForREU)
    {
        cout << "满足停机条件1：REQ-value≤varepsilon1 且 REU-value≤varepsilon2，分别为" << REQ << " " << REU << endl;
        return true;
    }
    else
    {
        //cout << "并不满足停机条件" << endl;
        return false;
    }
}

bool GenGroup::evolution()
{
    //Firstly set the objects in the parent generation
    int ParentsIndex[dimension];
    //cout << "Is generating the ParentIndex vector:" << endl;
    for (int i = 0; i < EliteNumber; i++)
    {
        ParentsIndex[i] = i;
        //cout << ParentsIndex[i]<<" ";
    }
    for (int i = EliteNumber; i < dimension; i++)
    {
        ParentsIndex[i] = rand() % (GroupScale - EliteNumber) + EliteNumber;
        //cout << ParentsIndex[i]<< " ";
    }
    for (int i = 0; i < ChildrenNumber; i++)
    {
        //cout <<endl<<"Is generating the " << i << "child" << endl;
        //cout << "Is generating the children gene of the " << i << "child" << endl;
        Gen* g = new Gen;
        double* p = new double[dimension];
        //cout << "Is generating the gene composition vector of the children gene of the " << i << "child" << endl;
        ChildrenSpaceGenerate(p);
        //cout << "It is the final q vector of this child:" << endl;
        for (int j = 0; j < x_matrix_len; j++)
        {
            g->q[j] = 0;
            for (int k = 0; k < dimension; k++)
            {
                g->q[j] += p[k]*this->Gens[ParentsIndex[k]].q[j];
            }
        }
        g->SmoothQ();
        g->GenUpdate();
        if (g->GenFitness<this->Gens[GroupScale-1].GenFitness)
        {
            Gens[GroupScale - 1].GenCopy(*g);
            //cout << "第" << i << "个子是优秀子，适应值为 "<<g->GenFitness<<" 可以进行更新" << endl;
            this->GensUpdate();
        }
        else
        {
            //cout << "第" << i << "个子不是优秀子，适应值为 " << g->GenFitness << "不进行更新" << endl;
        }
        delete[] p;
        delete g;
        if (CheckStop())
        {
            return true;
        }
    }

/*

    if (this->Gens[0].GenFitness == this->Gens[GroupScale - 1].GenFitness)
    {
        if (this->Gens[0]==this->Gens[GroupScale-1])
        {
            cout << "全局收敛了";
            return true;
        }
    }

*/
    return false;
}

bool GenGroup::variation()
{
    Gen* g_variation_group = new Gen[variation_num];
    for (int i = 0; i < variation_num; i++)
    {
        g_variation_group[i].GenCopy(Gens[i]);//将i改为0，突变更有针对性;但是会导致过快收敛
        double movement = (rand() % variation_gap_num_scale + 1) * alpha_gap;
        g_variation_group[i].q[(rand() % x_matrix_len)] += (rand() % 2 == 0) ? movement : (-1) * movement;
        //g_variation_group[i].SmoothQ();
        g_variation_group[i].GenUpdate();
    }
    for (int i = 0; i < variation_num; i++)
    {
        if (g_variation_group[i].GenFitness < Gens[GroupScale-1-i].GenFitness)
        {
            Gens[GroupScale-1-i].GenCopy(g_variation_group[i]);//发现gensupdate耗时很久，故作此优化
            //cout << "第" << i << "个变异子是优秀子，适应值为 "<< g_variation_group[i].GenFitness<<" 可以进行更新" << endl;

        }
        else
        {
            //cout << "第" << i << "个子不是优秀子，适应值为 " << g_variation_group[i].GenFitness << "不进行更新" << endl;
        }

    }
    this->GensUpdate();
    delete[] g_variation_group;
    if (CheckStop())
    {
        return true;
    }
    return false;
}
