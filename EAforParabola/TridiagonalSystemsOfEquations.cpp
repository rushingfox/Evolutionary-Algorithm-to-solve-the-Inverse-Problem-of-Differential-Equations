//׷��������ԽǷ�����
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
			cout << "���Խ�Ԫ����0�������޽⣡��" << endl;
			return b;
		}
		A[j + 1] = A[j + 1] / s;//һ���������ǹ�һ������Ԫ���Խ����µ�Ԫ��û��Ҫ���в���
		b[k] = b[k] / s;//��Ϊ����Ļش�û�õ���Щ����
		A[j + 3] = A[j + 3] - A[j + 2] * A[j + 1];
		b[k + 1] = b[k + 1] - A[j + 2] * b[k];
	}
	s = A[3 * t - 3];
	if (fabs(s) + 1.0 == 1.0)
	{
		cout << "��Ԫ�����Խ�Ԫ����0���޽⣡��" << endl;
		return b;
	}
	b[t - 1] = b[t - 1] / s;
	for (k = t - 2; k >= 0; k--)//�ش���ⷽ�̣����֮�ºܼ��ˣ��ۻ���Ͷ�û��Ҫ��
		b[k] = b[k] - b[k + 1] * A[3 * k + 1];
	return b;
}