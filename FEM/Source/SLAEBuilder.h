#pragma once

#include "FEMInfo.h"
#include "Math.h"
#include "Matrix.h"

#include <vector>
#include <algorithm>
#include <functional>

using namespace std;

class SLAEBuilder
{
public:
	SLAEBuilder(vector<Point>& points, vector<FiniteElement>& elements) :
		points(points),
		elements(elements)
	{
		local.resize(basisSize, vector<double>(basisSize));
		localb.resize(basisSize);

		BuildPatterns();
	}

	bool Build(Matrix& A, vector<double>& b)
	{
		if (lambda == nullptr || gamma == nullptr || f == nullptr)
			return false;

		for (auto& e : elements)
		{
			BuildLocalMatrix(e);
			BuildLocalB(e);
			AddLocalToGlobal(A, b, e);
		}

		return true;
	}

	void SetLambda(vector<double>* lambda) { this->lambda = lambda; }
	void SetGamma(vector<double>* gamma) { this->gamma = gamma; }
	void SetF(function<double(double, double)>* f) { this->f = f; }

private:
	vector<Point>& points;
	vector<FiniteElement>& elements;

	vector<double>* lambda = nullptr;
	vector<double>* gamma = nullptr;
	function<double(double, double)>* f = nullptr;

	vector<vector<double>> M;
	vector<vector<vector<LocalComp>>> gPattern;

	vector<vector<double>> local;
	vector<double> localb;

private:
	void BuildPatterns()
	{
		gPattern.resize(basisSize, vector<vector<LocalComp>>(basisSize));
		M.resize(basisSize, vector<double>(basisSize));

		for (int i = 0; i < basisSize; i++)
			for (int j = 0; j < basisSize; j++)
			{
				for (auto& x : ders[i])
					for (auto& y : ders[j])
					{
						int v1 = x.v1 + y.v1;
						int v2 = x.v2 + y.v2;
						int v3 = x.v3 + y.v3;

						double coeff = x.coeff * y.coeff * (double)Factorial(v1) * Factorial(v2) * Factorial(v3) / (double)Factorial(v1 + v2 + v3 + 2);
						gPattern[i][j].push_back(LocalComp(x.gradNo, y.gradNo, coeff));
					}

				double sum = 0;
				for (auto& x : basis[i])
					for (auto& y : basis[j])
					{
						int v1 = x.v1 + y.v1;
						int v2 = x.v2 + y.v2;
						int v3 = x.v3 + y.v3;

						sum += x.coeff * y.coeff * (double)Factorial(v1) * Factorial(v2) * Factorial(v3) / (double)Factorial(v1 + v2 + v3 + 2);
					}
				M[i][j] = sum;
			}
	}

	void BuildLocalMatrix(FiniteElement& e)
	{
		double D = abs(Det(e, points));
		vector<double> alpha = Alpha(e, points);
		vector<Grad> grads =
		{
			{ alpha[0], alpha[1] },
			{ alpha[2], alpha[3] },
			{ alpha[4], alpha[5] }
		};


		for (int i = 0; i < basisSize; i++)
			for (int j = 0; j < basisSize; j++)
			{
				double G = 0;
				for (auto comp : gPattern[i][j])
				{
					double scalGrad = grads[comp.grad1].a1 * grads[comp.grad2].a1 + grads[comp.grad1].a2 * grads[comp.grad2].a2;
					G += comp.coeff * scalGrad;
				}

				local[i][j] = (gamma->at(e.materialNo) * M[i][j] + lambda->at(e.materialNo) * G) * D;
			}
	}

	void BuildLocalB(FiniteElement& e)
	{
		vector<Point> coords = CalculateCoords(e, points);
		double D = abs(Det(e, points));

		for (int i = 0; i < basisSize; i++)
		{
			localb[i] = Triangle([=](double ksi, double etta)
				{
					double L1 = 1 - ksi - etta;
					double L2 = ksi;
					double L3 = etta;

					double res = 0.0;
					for (auto comp : basis[i])
						res += comp.coeff * pow(L1, comp.v1) * pow(L2, comp.v2) * pow(L3, comp.v3);

					double x = L1 * coords[0].x + L2 * coords[1].x + L3 * coords[2].x;
					double y = L1 * coords[0].y + L2 * coords[1].y + L3 * coords[2].y;

					res *= (*f)(x, y);

					return res;

				}, D);
		}
	}

	void AddLocalToGlobal(Matrix& A, vector<double>& b, FiniteElement& e)
	{
		for (int i = 0; i < basisSize; i++)
		{
			A.DI[e.verts[i]] += local[i][i];
			b[e.verts[i]] += localb[i];

			for (int j = 0; j < i; j++)
			{
				auto a = e.verts[i];
				auto b = e.verts[j];
				if (a < b) swap(a, b);

				auto begin = A.JA.begin() + A.IA[a];
				if (A.IA[a + 1] > A.IA[a])
				{
					auto end = A.JA.begin() + A.IA[a + 1] - 1;
					auto iter = lower_bound(begin, end, b);
					auto index = iter - A.JA.begin();
					A.AL[index] += local[i][j];
				}
			}
		}
	}
};