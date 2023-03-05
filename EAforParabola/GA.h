#ifndef GA_H
#define GA_H
#include <stdint.h>
#include <time.h>
#define _USE_MATH_DEFINES
#include <math.h>

//the boundary of searching
const double alpha_up = 10;
const double alpha_low = 0.1;
const double alpha_gap = 0.001;
const int alpha_max_num = (alpha_up - alpha_low) / alpha_gap;
//the number of the hat functions in the domain of definition-1(0-n)
const int n = 100;
//based n, give the length of some matrixes
const int x_matrix_len = n + 1;
const int b_and_u_matrix_len = n - 1;//b and u is the same length
const int A_matrix_len = 3 * n - 5;
//the max and min of x, the little gap in the hat function
const double x_min = 0;
const double x_max = 1;
const double x_gap = (double)((x_max - x_min) / n);

const int GroupScale = 200;
const int dimension = 10;
const int EliteNumber = 1;
const int ChildrenNumber = 1;
//controllment parameter for end
const double varepsilonForREQ = 0.01;
const double varepsilonForREQ_variation = 0.002;
const double varepsilonForREU = 0.05;
const double varepsilonForREU_variation = 0.008;

//regularization parameter
extern double beta;
//noise parameter
const double delta = 0.1;
const int noise_n = 100;//n等分
const double noise_gap = ((double)2/noise_n);//每段长度
const int noise_point_num = noise_n + 1;//可选点数量

//parameters for the children space generation
const double a_low = -0.5;
const double a_up = 1.5;
const double a_gap = 0.01;
const int a_n = ((a_up - a_low) / a_gap)+1;

extern double x_matrix[x_matrix_len];
//构造对称三对角方程的输入矩阵b
extern double b_matrix[b_and_u_matrix_len];
//构造对称三对角方程的输入矩阵A
extern double A_matrix[A_matrix_len];
extern double u_observe[b_and_u_matrix_len];

//T为时间区域的边界，也即我们用来反演的数据时间点
const double T = 1;

//para for variation
extern double variation_rate;
extern double variation_gap_num_scale_rate;
extern int variation_gap_num;
extern int variation_num;

class Gen
{

public:
    double q_value(double x);
    double q[x_matrix_len];
    double u[b_and_u_matrix_len];
    double GenFitness;
    double GenUpdate();
    void GenRandomGenerate();
    bool operator < (const Gen& s);
    void GenCopy(const Gen& g_origin);//不包含update操作，故建议先update（更新fitness）再使用此函数
    bool operator == (const Gen& g);
    void SmoothQ();
    void GenPrint();
};

class GenGroup
{
public:
    GenGroup();
    double REQ;
    double REU;
    Gen Gens[GroupScale];
    bool evolution();
    bool variation();//针对优秀个体进行细微变异
    void GensUpdate();
    void ChildrenSpaceGenerate(double * p);
    bool CheckStop();
};

void A_matrix_change(Gen* g);
double u(double x, double t);
double f(double x, double t);
double q(double x);
double fitness(Gen* g);
double EuclidDistance(double* m1, double* m2, int size);
#endif