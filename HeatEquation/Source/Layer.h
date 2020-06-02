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
		c = 1.0 / (t[1] - t[0]);
		cs.push_back(c);
	}

	void SetQ(vector<vector<double>*> q)
	{
		qs.clear();
		qs.push_back(q[0]);
	}
};

class ThreeLayer : public Layer
{
public:
	ThreeLayer() : Layer(2) {}

	void SetT(vector<double> t)
	{
		cs.clear();
		
		dt = t[2] - t[0];
		dt0 = t[2] - t[1];
		dt1 = t[1] - t[0];

		c = (dt + dt0) / (dt * dt0);

		cs.push_back(-dt0 / (dt * dt1));
		cs.push_back(dt / (dt1 * dt0));

	}

	void SetQ(vector<vector<double>*> q)
	{
		qs.clear();
		qs.push_back(q[0]);
		qs.push_back(q[1]);
	}

private:
	double dt, dt0, dt1;
};

class FourLayer : public Layer
{
public:
	FourLayer() : Layer(3) {}

	void SetT(vector<double> t)
	{
		cs.clear();

		dt03 = t[3] - t[0];
		dt02 = t[3] - t[1];
		dt01 = t[3] - t[2];

		dt13 = t[2] - t[0];
		dt12 = t[2] - t[1];

		dt23 = t[1] - t[0];

		c = (dt02 * dt01 + dt03 * dt01 + dt03 * dt02) / (dt03 * dt02 * dt01);

		cs.push_back(dt02 * dt01 / (dt23 * dt13 * dt03));
		cs.push_back(-dt03 * dt01 / (dt23 * dt12 * dt02));
		cs.push_back(dt03 * dt02 / (dt13 * dt12 * dt01));
	}

	void SetQ(vector<vector<double>*> q)
	{
		qs.clear();
		qs.push_back(q[0]);
		qs.push_back(q[1]);
		qs.push_back(q[2]);
	}

private:
	double dt01, dt02, dt03, dt12, dt13, dt23;
};


