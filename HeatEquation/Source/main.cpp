#include "FEMInfo.h"
#include "Math.h"
#include "Matrix.h"

#include "MeshBuilder.h"
#include "Interval.h"
#include "TimeMeshBuilder.h"
#include "PortraitBuilder.h"
#include "SLAEBuilder.h"

#include "Solver.h"

#include "Layer.h"

#include <string>
#include <fstream>

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

void InputTime(string timeFile, vector<Interval>& time_intervals)
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
//void InputBound(string bounds1File, string bounds2File)
//{
//	ifstream in(bounds1File);
//	in >> edgeCount1;
//	bound1.resize(edgeCount1);
//	for (int i = 0; i < edgeCount1; i++)
//	{
//		in >> bound1[i].v1 >> bound1[i].v4 >> bound1[i].valueNo;
//		auto a = bound1[i].v1;
//		auto b = bound1[i].v4;
//		bool f = a > b;
//		if (f) swap(a, b);
//
//		bound1[i].v2 = edgeMatrix[a][b] + (f ? 0 : 1);
//		bound1[i].v3 = edgeMatrix[a][b] + (f ? 1 : 0);
//	}
//	in.close();
//
//	in.open(bounds2File);
//	in >> edgeCount2;
//	bound2.resize(edgeCount2);
//	for (int i = 0; i < edgeCount2; i++)
//	{
//		in >> bound2[i].v1 >> bound2[i].v4 >> bound2[i].valueNo;
//		auto a = bound2[i].v1;
//		auto b = bound2[i].v4;
//		bool f = a > b;
//		if (f) swap(a, b);
//
//		bound2[i].v2 = edgeMatrix[a][b] + (f ? 0 : 1);
//		bound2[i].v3 = edgeMatrix[a][b] + (f ? 1 : 0);
//	}
//	in.close();
//}


int main()
{
	int nodeCount;
	vector<Point> points;
	vector<FiniteElement> elements;

	InputGrid(nodeCount, "points.txt", points, "triangles.txt", elements);
	MeshBuilder mesh(nodeCount);
	mesh.Build(elements);
	nodeCount = mesh.nodeCount;

	vector<Interval> time_intervals;
	InputTime("time.txt", time_intervals);
	TimeMeshBuilder time(time_intervals);

	const vector<double> t = time.getMesh();

	Matrix A;
	PortraitBuilder portrait(nodeCount, elements);
	portrait.Build(A);

	vector<double> b(nodeCount);

	SLAEBuilder builder(points, elements);

	vector<double> q0(nodeCount), q1(nodeCount), q2(nodeCount);
	vector<vector<double>*> Qs;

	for (int i = 0; i < A.N; i++)
		q0.push_back(u(points[i].x, points[i].y, t[0]));

	Qs.push_back(&q0);
	Layer* layer = new TwoLayer();

	layer->SetQ({ Qs[0] });
	layer->SetT({ t[0], t[1] });

	builder.SetLayer(layer);
	builder.Build(A, b);

	// Учет краевых

	Solvers::BCG(A, q1, b);
	
	A.Clear();

	Qs.push_back(&q1);

	Layer* layer = new ThreeLayer();

	layer->SetQ({ Qs[0], Qs[1] });
	layer->SetT({ t[0], t[1], t[2] });

	builder.SetLayer(layer);
	builder.Build(A, b);

	// Учет краевых

	Solvers::BCG(A, q2, b);

	A.Clear();

	Qs.push_back(&q2);

	Layer* layer = new FourLayer();
	for (int i = 3; i < t.size(); i++)
	{
		layer->SetQ({ Qs[i - 2], Qs[i - 1], Qs[i] });
		layer->SetT({ t[i - 3], t[i - 2], t[i - 1], t[i] });

		builder.SetLayer(layer);
		builder.Build(A, b);

		vector<double> q(nodeCount);

		// Учет краевых

		Solvers::BCG(A, q, b);

		A.Clear();

		Qs.push_back(&q);
	}

	system("pause");
	return 0;
}