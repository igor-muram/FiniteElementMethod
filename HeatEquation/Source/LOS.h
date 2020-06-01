#pragma once

#include <vector>

using namespace std;

//void Compute()
//{
//	vector<double> r, z, Az;
//	int maxiter = 100000;
//	double eps = 1.E-30;
//	r.resize(nodeCount);
//	z.resize(nodeCount);
//	Az.resize(nodeCount);
//
//	// r = A*x
//	Multiply(q, r);
//
//	// r = f - A*x 
//	for (int i = 0; i < nodeCount; i++)
//	{
//		r[i] = globalB[i] - r[i];
//		z[i] = r[i];
//	}
//
//	double normB = 1.0;
//	double scal_rr = Scal(r, r);
//	double diff = sqrt(fabs(scal_rr)) / normB;
//	int iter = 0;
//	for (; iter < maxiter && diff >= eps; iter++)
//	{
//		Multiply(z, Az);     // Az = A*z
//		double a = scal_rr / Scal(Az, z);
//		for (int i = 0; i < nodeCount; i++)
//		{
//			q[i] += a * z[i];
//			r[i] -= a * Az[i];
//		}
//
//		double b = 1. / scal_rr;
//		scal_rr = Scal(r, r);
//		b *= scal_rr;
//		for (int i = 0; i < nodeCount; i++)
//			z[i] = r[i] + b * z[i];
//
//		diff = sqrt(fabs(Scal(r, r))) / normB;
//	}
//}
//
//void Multiply(vector<double>& x, vector<double>& res)
//{
//	for (int i = 0; i < nodeCount; i++)
//	{
//		res[i] = globalMatrix.DI[i] * x[i];
//		int i0 = globalMatrix.IA[i];
//		int i1 = globalMatrix.IA[i + 1];
//		for (int k = i0; k < i1; k++)
//		{
//			int j = globalMatrix.JA[k];
//			res[i] += globalMatrix.AL[k] * x[j];
//			res[j] += globalMatrix.AL[k] * x[i];
//		}
//	}
//}
//
//double FEM::Scal(vector<double>& x, vector<double>& y)
//{
//	double scal = 0.0;
//	for (int i = 0; i < nodeCount; i++)
//		scal += x[i] * y[i];
//
//	return scal;
//}