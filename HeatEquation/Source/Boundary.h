#pragma once

#include "FEMInfo.h"

struct Edge
{
	int v1, v2, v3, v4;
	int valueNo;
};

void Boundary1(Matrix& A, vector<double>& b, vector<Edge>& bound, map<int, Point>& pointsMap, double t)
{
	for (auto edge : bound)
	{
		A.DI[edge.v1] = 1.0E+50;
		b[edge.v1] = 1.0E+50 * uValue[edge.valueNo](pointsMap[edge.v1].x, pointsMap[edge.v1].y, t);
		A.DI[edge.v2] = 1.0E+50;
		b[edge.v2] = 1.0E+50 * uValue[edge.valueNo](pointsMap[edge.v2].x, pointsMap[edge.v2].y, t);
		A.DI[edge.v3] = 1.0E+50;
		b[edge.v3] = 1.0E+50 * uValue[edge.valueNo](pointsMap[edge.v3].x, pointsMap[edge.v3].y, t);
		A.DI[edge.v4] = 1.0E+50;
		b[edge.v4] = 1.0E+50 * uValue[edge.valueNo](pointsMap[edge.v4].x, pointsMap[edge.v4].y, t);
	}
}

//void Boundary2(Matrix& A, vector<double>& b, vector<Edge>& bound, map<int, Point>& pointsMap, double t)
//{
//	vector<double> conditions(4, 0.0);
//	vector<vector<double>> M =
//	{
//		{ 8.0 / 105, 19.0 / 1680, 33.0 / 560, -3.0 / 140 },
//		{ 19.0 / 1680, 8.0 / 105, -3.0 / 140, 33.0 / 560 },
//		{  33.0 / 560, -3.0 / 140, 27.0 / 70, -27.0 / 560 },
//		{  -3.0 / 140, 33.0 / 560, -27.0 / 560, 27.0 / 70 }
//	};
//
//	for (int edge = 0; edge < bound.size(); edge++)
//	{
//		vector<int> v = { bound[edge].v1, bound[edge].v2, bound[edge].v3, bound[edge].v4 };
//		int thetaNo = bound[edge].valueNo;
//		
//		vector<double> theta(4);
//		for (int i = 0; i < 4; i++)
//			theta[i] = thetaValue[thetaNo](pointsMap[v[i]].x, pointsMap[v[i]].y, t);
//
//		double h = Distance(pointsMap[v[0]], pointsMap[v[3]]);
//
//		fill(conditions.begin(), conditions.end(), 0.0);
//		for (int i = 0; i < 4; i++)
//			for (int j = 0; j < 4; j++)
//				conditions[i] += theta[j] * M[i][j];
//
//		b[v[0]] += conditions[0];
//		b[v[1]] += conditions[1];
//		b[v[2]] += conditions[2];
//		b[v[3]] += conditions[3];
//	}
//}

void Boundary2(double lambda, Matrix& A, vector<double>& b, vector<Edge>& bound, map<int, Point>& pointsMap, double t)
{
	vector<double> conditions(4, 0.0);


	for (int edge = 0; edge < bound.size(); edge++)
	{
		vector<int> v = { bound[edge].v1, bound[edge].v2, bound[edge].v3, bound[edge].v4 };
		int thetaNo = bound[edge].valueNo;

		function<double(double)> theta = thetaValue[thetaNo];
		double h = Distance(pointsMap[v[0]], pointsMap[v[3]]);

		vector<function<double(double)>> basis =
		{
			[h](double x) { return 0.5 * x * (3 * x - 1) * (3 * x - 2); },
			[h](double x) { return 4.5 * (1 - x) * x * (3 * x - 1); },
			[h](double x) { return 4.5 * (1 - x) * x * (3 * (1 - x) - 1); },
			[h](double x) { return 0.5 * (1 - x) * (3 * (1 - x) - 1) * (3 * (1 - x) - 2); }
		};

		for (int i = 0; i < 4; i++)
		{
			function<double(double)> psi = basis[i];
			conditions[i] = h * lambda * NewtonCotes(0.0, 1, [theta, psi](double x) { return theta(x) * psi(x); });
		}

		b[v[0]] += conditions[0];
		b[v[1]] += conditions[1];
		b[v[2]] += conditions[2];
		b[v[3]] += conditions[3];
	}
}