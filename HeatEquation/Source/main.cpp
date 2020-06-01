#include "FEMInfo.h"
#include "Math.h"
#include "Matrix.h"

#include "MeshBuilder.h"
#include "PortraitBuilder.h"
#include "SLAEBuilder.h"

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

	Matrix A;
	PortraitBuilder portrait(nodeCount, elements);
	portrait.Build(A);

	vector<double> b(nodeCount);

	SLAEBuilder builder(points, elements);

	vector<double> q0(nodeCount);
	vector<vector<double>*> Qs;
	Qs.push_back(&q0);
	Layer* layer = new TwoLayer();

	layer->SetQ({ &q0 });
	layer->SetT({ 1.0 });

	builder.SetLayer(layer);
	builder.Build(A, b);
	
	A.Clear();

	system("pause");
	return 0;
}