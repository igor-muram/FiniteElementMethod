#pragma once

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