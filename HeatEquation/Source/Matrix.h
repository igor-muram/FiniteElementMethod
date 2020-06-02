#pragma once

#include <vector>
#include <stdexcept>

using namespace std;

struct Matrix
{
	int N;
	vector<int> IA, JA;
	vector<double> DI, AL;

	double& operator()(int i, int j)
	{
		if (i >= N || j >= N || i < 0 || j < 0)
			throw std::out_of_range("Bad index for matrix");

		if (i == j)
			return DI[i];

		if (j > i)
		{
			int width = IA[j + 1] - IA[j];

			int begin = j - width;
			int end = j;
			for (int k = begin; k < end; k++)
			{
				if (k == i)
				{
					return AL[IA[j] + k - begin];
				}
			}

		}
		else
		{
			int width = IA[i + 1] - IA[i];

			int begin = i - width;
			int end = i;
			for (int k = begin; k < end; k++)
			{
				if (k == j)
				{
					return AL[IA[i] + k - begin];
				}
			}
		}

		throw std::out_of_range("Bad index for matrix");
	}

	void Clear()
	{
		fill(DI.begin(), DI.end(), 0.0);
		fill(AL.begin(), AL.end(), 0.0);
	}
};
