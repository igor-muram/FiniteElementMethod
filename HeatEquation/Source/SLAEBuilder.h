#pragma once

#include "FEMInfo.h"
#include "Math.h"
#include "Matrix.h"

#include "Layer.h"

#include <vector>
#include <algorithm>
using namespace std;

class SLAEBuilder
{
public:
	SLAEBuilder(vector<Point>& points, vector<FiniteElement>& elements, int nodeCount) : 
		nodeCount(nodeCount),
		points(points), 
		elements(elements), 
		layer(layer)
	{
		local.resize(basisSize, vector<double>(basisSize));
		localb.resize(basisSize);

		BuildPatterns();
	}

	void Build(Matrix& A, vector<double>& b)
	{
			for (auto& e : elements)
			{
				BuildLocalMatrix(e);
				BuildLocalB(e);
				AddLocalToGlobal(A, b, e);
			}
	}

	void Boundary(Matrix& A, vector<double>& b, double u1, double u2)
	{
		A(0, 0) = 1.0e+50;
		b[0] = u1 * 1.0e+50;

		A(nodeCount - 1, nodeCount - 1) = 1.0e+50;
		b[nodeCount - 1] = u2 * 1.0e+50;
	}

	void SetLayer(Layer* layer)
	{
		this->layer = layer;
	}

private:
	vector<Point>& points;
	vector<FiniteElement>& elements;
	int nodeCount;

	Layer* layer;

	vector<vector<double>> local;
	vector<double> localb;

	vector<vector<double>> M;
	vector<vector<vector<LocalComp>>> gPattern;

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

		double C = layer->c;

		for (int i = 0; i < basisSize; i++)
			for (int j = 0; j < basisSize; j++)
			{
				double G = 0;
				for (auto comp : gPattern[i][j])
				{
					double scalGrad = grads[comp.grad1].a1 * grads[comp.grad2].a1 + grads[comp.grad1].a2 * grads[comp.grad2].a2;
					G += comp.coeff * scalGrad;
				}
				local[i][j] = (C * gamma[e.materialNo] * M[i][j] + lambda[e.materialNo] * G) * D;
			}
	}

	void BuildLocalB(FiniteElement& e)
	{
		vector<Point> coords = CalculateCoords(e, points);
		double D = abs(Det(e, points));

		const vector<vector<double>*> Qs = layer->qs;
		const vector<double> Cs = layer->cs;
		int size = layer->size;

		vector<double> temp(basisSize);
		for (int i = 0; i < basisSize; i++)
			temp[i] = f[e.materialNo](coords[i].x, coords[i].y, 0.0);

		for (int i = 0; i < basisSize; i++)
		{
			localb[i] = temp[0] * M[i][0];

			for (int j = 1; j < basisSize; j++)
				localb[i] += temp[j] * M[i][j];

			for (int i = 0; i < size; i++)
				for (int j = 0; j < basisSize; j++)
					localb[i] += Cs[i] * (*Qs[i])[e.verts[j]] * M[i][j];

			localb[i] *= D;
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