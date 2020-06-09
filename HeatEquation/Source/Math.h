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

	for (double x = 0.0; x < b; x += h)
		result += 32 * f(x + h / 4) + 12 * f(x + h / 2) + 32 * f(x + 3 * h / 4);

	result = result * 0.5 * h / 45;

	return result;
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