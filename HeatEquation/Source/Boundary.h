#pragma once

#include "FEMInfo.h"

struct Edge
{
	int v1, v2, v3, v4;
	int valueNo;
};

void Boundary1(Matrix& A, vector<double>& b, vector<Edge>& bound, map<int, Point>& pointsMap, double t)
{
	vector<vector<double>> M =
	{
		{  8.0 / 105, 33.0 / 560, -3.0 / 140, 19.0 / 1680 },
		{ 33.0 / 560, 27.0 / 70, -27.0 / 560, -3.0 / 140 },
		{ -3.0 / 140, -27.0 / 560, 27.0 / 70, 33.0 / 560 },
		{ 19.0 / 1680, -3.0 / 140, 33.0 / 560, 8.0 / 105 }
	};

	vector<function<double(double)>> basis =
	{
		[](double x) { return 0.5 * (1 - x) * (3 * (1 - x) - 1) * (3 * (1 - x) - 2); },
		[](double x) { return 4.5 * (1 - x) * x * (3 * (1 - x) - 1); },
		[](double x) { return 4.5 * (1 - x) * x * (3 * x - 1); },
		[](double x) { return 0.5 * x * (3 * x - 1) * (3 * x - 2); }
	};

	for (auto edge : bound)
	{
		function<double(double, double)> Ug = uValue[edge.valueNo];
		vector<double> f(4);

		double x0 = pointsMap[edge.v1].x;
		double x1 = pointsMap[edge.v4].x;
		double y0 = pointsMap[edge.v1].y;
		double y1 = pointsMap[edge.v4].y;

		for (int i = 0; i < 4; i++)
			f[i] = NewtonCotes(0.0, 1.0, [&](double t) { return Ug(x0 + t * (x1 - x0), y0 + t * (y1 - y0)) * basis[i](t); });

		vector<double> q(4);
		Gauss(M, q, f);

		A.DI[edge.v1] = 1.0E+50;
		b[edge.v1] = 1.0E+50 * q[0];
		A.DI[edge.v2] = 1.0E+50;
		b[edge.v2] = 1.0E+50 * q[1];
		A.DI[edge.v3] = 1.0E+50;
		b[edge.v3] = 1.0E+50 * q[2];
		A.DI[edge.v4] = 1.0E+50;
		b[edge.v4] = 1.0E+50 * q[3];
	}
}

void Boundary2(Matrix& A, vector<double>& b, vector<Edge>& bound, map<int, Point>& pointsMap, double lambda, double t)
{
	vector<function<double(double)>> basis =
	{
		[](double x) { return 0.5 * (1 - x) * (3 * (1 - x) - 1) * (3 * (1 - x) - 2); },
		[](double x) { return 4.5 * (1 - x) * x * (3 * (1 - x) - 1); },
		[](double x) { return 4.5 * (1 - x) * x * (3 * x - 1); },
		[](double x) { return 0.5 * x * (3 * x - 1) * (3 * x - 2); }
	};

	for (int edge = 0; edge < bound.size(); edge++)
	{
		vector<int> v = { bound[edge].v1, bound[edge].v2, bound[edge].v3, bound[edge].v4 };
		int thetaNo = bound[edge].valueNo;

		function<double(double)> theta = thetaValue[thetaNo];
		double h = Distance(pointsMap[v[0]], pointsMap[v[3]]);

		double x1 = pointsMap[v[0]].x;
		double x2 = pointsMap[v[3]].x;

		for (int i = 0; i < 4; i++)
		{
			function<double(double)> psi = basis[i];
			b[v[i]] += h * NewtonCotes(0.0, 1.0, [&](double x) { return theta(x1 + (x2 - x1) * x) * psi(x); });
		}

	}
}