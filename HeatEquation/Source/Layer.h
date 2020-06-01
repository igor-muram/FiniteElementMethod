#pragma once

#include "SLAEBuilder.h"

#include <vector>

using namespace std;

class Layer
{
	friend class SLAEBuilder;

public:
	Layer(int size) : size(size)
	{};

	virtual void SetT(vector<double> t) = 0;
	virtual void SetQ(vector<vector<double>*> q) = 0;

protected:
	vector<vector<double>*> qs;
	vector<double> cs;
	double c;

	int size;
};

class TwoLayer : public Layer
{
public:
	TwoLayer() : Layer(1) {}

	void SetT(vector<double> t)
 	{
		cs.clear();
		cs.push_back(t[0]);

		c = t[0];
	}

	void SetQ(vector<vector<double>*> q)
	{
		qs.clear();
		qs.push_back(q[0]);
	}
};

