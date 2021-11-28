#include "FEMInfo.h"
#include "Math.h"
#include "Matrix.h"
#include "Input.h"
#include "MeshBuilder.h"
#include "Interval.h"
#include "PortraitBuilder.h"
#include "SLAEBuilder.h"
#include "Boundary.h"
#include "Solver.h"

#include <string>
#include <fstream>

using namespace std;

bool PointInsideTriangle(Point t1, Point t2, Point t3, Point p)
{
	double crossProduct1 = (t1.x - p.x) * (t2.y - t1.y) - (t2.x - t1.x) * (t1.y - p.y);
	double crossProduct2 = (t2.x - p.x) * (t3.y - t2.y) - (t3.x - t2.x) * (t2.y - p.y);
	double crossProduct3 = (t3.x - p.x) * (t1.y - t3.y) - (t1.x - t3.x) * (t3.y - p.y);

	if (crossProduct1 >= 0.0 && crossProduct2 >= 0.0 && crossProduct3 >= 0.0)
		return true;

	if (crossProduct1 <= 0.0 && crossProduct1 <= 0.0 && crossProduct3 <= 0.0)
		return true;

	return false;
}

int main()
{
	int nodeCount;
	vector<Point> points;
	vector<FiniteElement> elements;

	InputGrid(nodeCount, "points.txt", points, "triangles.txt", elements);
	MeshBuilder mesh(nodeCount);
	mesh.Build(elements);
	nodeCount = mesh.nodeCount;

	vector<Edge> bound1, bound2;
	InputBound(bound1, "boundary1.txt");
	mesh.BuildBoundary(bound1);
	InputBound(bound2, "boundary2.txt");
	mesh.BuildBoundary(bound2);

	Matrix A;
	PortraitBuilder portrait(nodeCount, elements);
	portrait.Build(A);

	SLAEBuilder builder(points, elements);
	vector<double> b(nodeCount);

	map<int, Point> pointsMap;
	CreatePointMap(elements, points, pointsMap);

	vector<double> lambda = { 1 };
	vector<double> gamma = { 0 };
	function<double(double, double)> f = [](double x, double y) { return -12 * x * x; };
	builder.SetLambda(&lambda);
	builder.SetGamma(&gamma);
	builder.SetF(&f);

	vector<double> q(nodeCount);

	builder.Build(A, b);
	Boundary1(A, b, bound1, pointsMap);

	Solvers::BCG(A, q, b);

	Point p = Point(1.0 / 4, 1.0 / 10);
	FiniteElement found;

	for (int i = 0; i < elements.size(); i++)
	{
		FiniteElement e = elements[i];
		Point t2 = pointsMap[e.verts[0]];
		Point t1 = pointsMap[e.verts[1]];
		Point t3 = pointsMap[e.verts[2]];

		if (PointInsideTriangle(t1, t2, t3, p)) {
			found = e;
			break;
		}
	}

	vector<double> L = Ls(p, pointsMap[found.verts[0]], pointsMap[found.verts[1]], pointsMap[found.verts[2]]);

	double res = 0.0;
	for (int i = 0; i < basisSize; i++)
	{
		res += basisValue(i, L) * q[found.verts[i]];
	}

	for (int i = 0; i < nodeCount; i++)
	{
		cout << q[i] << endl;
	}

	cout << endl << endl << endl;

	cout << "points" << endl;

	for (int i = 0; i < pointsMap.size(); i++)
	{
		cout << pointsMap[i].x << ' ' << pointsMap[i].y << endl;
	}


	cout << res << endl;

	system("pause");
	return 0;
}