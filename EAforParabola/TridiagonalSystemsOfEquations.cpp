//追赶求解三对角方程组
#include  <iostream>
#include  <fstream>
#include  <cmath>
using namespace std;

double* trde(int t, double* A, double* b)
{
	int k, j;
	double s;
	for (k = 0; k < t - 1; k++)
	{
		j = 3 * k; s = A[j];
		if (fabs(3) + 1.0 == 1.0)
		{
			cout << "主对角元素是0，方程无解！！" << endl;
			return b;
		}
		A[j + 1] = A[j + 1] / s;//一下这四行是归一化和消元，对角线下的元素没必要进行操作
		b[k] = b[k] / s;//因为后面的回代没用到这些数据
		A[j + 3] = A[j + 3] - A[j + 2] * A[j + 1];
		b[k + 1] = b[k + 1] - A[j + 2] * b[k];
	}
	s = A[3 * t - 3];
	if (fabs(s) + 1.0 == 1.0)
	{
		cout << "消元后，主对角元素是0，无解！！" << endl;
		return b;
	}
	b[t - 1] = b[t - 1] / s;
	for (k = t - 2; k >= 0; k--)//回代求解方程，相比之下很简单了，累积求和都没必要的
		b[k] = b[k] - b[k + 1] * A[3 * k + 1];
	return b;
}