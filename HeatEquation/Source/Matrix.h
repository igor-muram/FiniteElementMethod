#pragma once

#include <vector>
#include <stdexcept>

using namespace std;

struct Matrix
{
	int N = 0;
	vector<int> IA, JA;
	vector<double> DI, AL;

	void Clear()
	{
		fill(DI.begin(), DI.end(), 0.0);
		fill(AL.begin(), AL.end(), 0.0);
	}
};
