#pragma once

#include <fstream>
#include <string>
#include <vector>

#include "Boundary.h"
#include "Interval.h"

using namespace std;

void InputGrid(
	int& nodeCount,
	string pointsFile,
	vector<Point>& points,
	string trianglesFile,
	vector<FiniteElement>& elements
)
{
	ifstream in(pointsFile);
	in >> nodeCount;
	points.resize(nodeCount);

	for (int i = 0; i < nodeCount; i++)
		in >> points[i].x >> points[i].y;

	in.close();

	in.open(trianglesFile);
	int triangleCount;
	in >> triangleCount;
	elements.resize(triangleCount);

	for (int i = 0; i < triangleCount; i++)
	{
		elements[i].verts.resize(basisSize);
		in >> elements[i].verts[0]
			>> elements[i].verts[1]
			>> elements[i].verts[2]
			>> elements[i].materialNo;
	}

	in.close();
}

void InputTime(vector<Interval>& time_intervals, string timeFile)
{
	ifstream in(timeFile);

	int time_count;
	in >> time_count;

	for (int i = 0; i < time_count; i++)
	{
		Interval interval;
		in >> interval;
		time_intervals.push_back(interval);
	}

	in.close();
}

void InputBound(vector<Edge>& bound, string boundsFile)
{
	int edgeCount;
	ifstream in(boundsFile);
	in >> edgeCount;
	bound.resize(edgeCount);

	for (int i = 0; i < edgeCount; i++)
		in >> bound[i].v1 >> bound[i].v4 >> bound[i].valueNo;

	in.close();
}