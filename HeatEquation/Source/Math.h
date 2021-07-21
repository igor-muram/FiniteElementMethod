#pragma once

#include <vector>
#include <map>
#include <set>
#include "FEMInfo.h"

using namespace std;

struct Point
{
	Point(double x, double y) : x(x), y(y) {}
	Point() : x(0.0), y(0.0) {}
	double x, y;

	Point operator+(Point& a)
	{
		return Point(a.x + x, a.y + y);
	}

	Point operator*(double constant)
	{
		return Point(x * constant, y * constant);
	}

	Point operator/(double constant)
	{
		return Point(x / constant, y / constant);
	}
};

struct Grad
{
	double a1, a2;
};

int Factorial(int n)
{
	int res = 1;
	for (int i = 2; i <= n; i++)
		res *= i;
	return res;
}

double Det(FiniteElement& e, vector<Point>& points)
{
	double x1 = points[e.verts[0]].x;
	double y1 = points[e.verts[0]].y;
	double x2 = points[e.verts[1]].x;
	double y2 = points[e.verts[1]].y;
	double x3 = points[e.verts[2]].x;
	double y3 = points[e.verts[2]].y;

	return (x2 - x1) * (y3 - y1) - (x3 - x1) * (y2 - y1);
}

double Det(Point& p1, Point& p2, Point& p3)
{
	double x1 = p1.x;
	double y1 = p1.y;
	double x2 = p2.x;
	double y2 = p2.y;
	double x3 = p3.x;
	double y3 = p3.y;

	return (x2 - x1) * (y3 - y1) - (x3 - x1) * (y2 - y1);
}

vector<double> Alpha(FiniteElement& e, vector<Point>& points)
{
	vector<double> alpha(6);
	double x1 = points[e.verts[0]].x;
	double y1 = points[e.verts[0]].y;
	double x2 = points[e.verts[1]].x;
	double y2 = points[e.verts[1]].y;
	double x3 = points[e.verts[2]].x;
	double y3 = points[e.verts[2]].y;

	double D = Det(e, points);

	alpha[0] = (y2 - y3) / D;
	alpha[1] = (x3 - x2) / D;
	alpha[2] = (y3 - y1) / D;
	alpha[3] = (x1 - x3) / D;
	alpha[4] = (y1 - y2) / D;
	alpha[5] = (x2 - x1) / D;

	return alpha;
}

vector<Point> CalculateCoords(FiniteElement& e, vector<Point>& points)
{
	Point t1 = points[e.verts[0]];
	Point t2 = points[e.verts[1]];
	Point t3 = points[e.verts[2]];

	vector<Point> coords;

	coords.push_back(t1);
	coords.push_back(t2);
	coords.push_back(t3);
	coords.push_back(((t2 / 2.0) + t1) * 2.0 / 3.0);
	coords.push_back(((t2 * 2.0) + t1) / 3.0);
	coords.push_back(((t3 / 2.0) + t1) * 2.0 / 3.0);
	coords.push_back(((t3 * 2.0) + t1) / 3.0);
	coords.push_back(((t3 / 2.0) + t2) * 2.0 / 3.0);
	coords.push_back(((t3 * 2.0) + t2) / 3.0);
	coords.push_back((t1 + t2 + t3) / 3.0);

	return coords;
}

double Distance(Point a, Point b)
{
	return sqrt((a.x - b.x) * (a.x - b.x) + (a.y - b.y) * (a.y - b.y));
}

void CreatePointMap(vector<FiniteElement>& elements, vector<Point>& points, map<int, Point>& pointsMap)
{
	for (auto e : elements)
	{
		Point t1 = points[e.verts[0]];
		Point t2 = points[e.verts[1]];
		Point t3 = points[e.verts[2]];

		pointsMap.try_emplace(e.verts[0], t1);
		pointsMap.try_emplace(e.verts[1], t2);
		pointsMap.try_emplace(e.verts[2], t3);
		pointsMap.try_emplace(e.verts[3], ((t2 / 2.0) + t1) * 2.0 / 3.0);
		pointsMap.try_emplace(e.verts[4], ((t2 * 2.0) + t1) / 3.0);
		pointsMap.try_emplace(e.verts[5], ((t3 / 2.0) + t1) * 2.0 / 3.0);
		pointsMap.try_emplace(e.verts[6], ((t3 * 2.0) + t1) / 3.0);
		pointsMap.try_emplace(e.verts[7], ((t3 / 2.0) + t2) * 2.0 / 3.0);
		pointsMap.try_emplace(e.verts[8], ((t3 * 2.0) + t2) / 3.0);
		pointsMap.try_emplace(e.verts[9], (t1 + t2 + t3) / 3.0);
	}
}

double NewtonCotes(double a, double b, function<double(double)> f)
{
	double h = (b - a) / 1000;
	double result = f(a);

	for (double x = h; x < b; x += h)
		result += 2 * f(x);

	result += f(b);
	result *= 7;

	for (double x = a; x < b; x += h)
		result += 32 * f(x + h / 4) + 12 * f(x + h / 2) + 32 * f(x + 3 * h / 4);

	result = result * 0.5 * h / 45;

	return result;
}

double Triangle(function<double(double, double)> f, double D)
{
	const int n = 21;
	vector<double> p1 = { 0.0451890097844, 0.0451890097844, 0.9096219804312, 0.7475124727339, 0.2220631655373, 0.7475124727339,
	0.2220631655373, 0.0304243617288, 0.0304243617288, 0.1369912012649, 0.6447187277637, 0.1369912012649, 0.2182900709714,
	0.2182900709714, 0.6447187277637, 0.0369603304334, 0.4815198347833, 0.4815198347833, 0.4036039798179, 0.4036039798179,
	0.1927920403641 };
	vector<double> p2 = { 0.0451890097844, 0.9096219804312, 0.0451890097844, 0.0304243617288, 0.0304243617288, 0.2220631655373,
	0.7475124727339, 0.7475124727339, 0.2220631655373, 0.2182900709714, 0.2182900709714, 0.6447187277637, 0.6447187277637,
	0.1369912012649, 0.1369912012649, 0.4815198347833, 0.0369603304334, 0.4815198347833, 0.1927920403641, 0.4036039798179,
	0.4036039798179 };
	vector<double> w = { 0.0519871420646, 0.0519871420646, 0.0519871420646, 0.0707034101784, 0.0707034101784, 0.0707034101784,
	0.0707034101784, 0.0707034101784, 0.0707034101784, 0.0909390760952, 0.0909390760952, 0.0909390760952, 0.0909390760952,
	0.0909390760952, 0.0909390760952, 0.1032344051380, 0.1032344051380, 0.1032344051380, 0.1881601469167, 0.1881601469167,
	0.1881601469167 };

	double sum = 0.0;

	for (int i = 0; i < n; i++)
		sum += f(p1[i], p2[i]) * w[i];

	return 0.25 * sum * D;
}

void Gauss(vector<vector<double>> A, vector<double>& x, vector<double> b)
{
	int N = A.size();
	for (int k = 0; k < N - 1; k++)
	{
		// Поиск ведущего элемента
		double max = abs(A[k][k]);
		int m = k;
		for (int i = k + 1; i < N; i++)
			if (abs(A[k][i]) > max)
			{
				max = abs(A[k][i]);
				m = i;
			}

		// Обмен местами b[m] и b[k]
		std::swap(b[m], b[k]);
		// Обмен местами k-ого и m-ого столбцов
		for (int j = k; j < N; j++)
			std::swap(A[k][j], A[m][j]);

		// Обнуление k-ого столбца
		for (int i = k + 1; i < N; i++)
		{
			double t = A[i][k] / A[k][k];
			b[i] -= t * b[k];
			for (int j = k + 1; j < N; j++)
				A[i][j] -= t * A[k][j];
		}
	}

	// Вычисление вектора x
	x[N - 1] = b[N - 1] / A[N - 1][N - 1];
	for (int k = N - 2; k >= 0; k--)
	{
		double sum = 0;
		for (int j = k + 1; j < N; j++)
			sum += A[k][j] * x[j];

		x[k] = (b[k] - sum) / A[k][k];
	}
}

vector<double> Ls(Point p, Point t1, Point t2, Point t3)
{
	vector<double> L(3);
	double D = abs(Det(t1, t2, t3));
	double D1 = abs(Det(p, t2, t3));
	double D2 = abs(Det(t1, p, t3));
	double D3 = abs(Det(t1, t2, p));

	L[0] = D1 / D;
	L[1] = D2 / D;
	L[2] = D3 / D;
	return L;
}

double basisValue(int i, vector<double> L)
{
	double res = 0.0;
	for (auto comp : basis[i])
		res += comp.coeff * pow(L[0], comp.v1) * pow(L[1], comp.v2) * pow(L[2], comp.v3);

	return res;
}