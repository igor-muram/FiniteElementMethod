#include "FEMInfo.h"
#include "Math.h"
#include "Matrix.h"
#include "Input.h"
#include "MeshBuilder.h"
#include "Interval.h"
#include "TimeMeshBuilder.h"
#include "PortraitBuilder.h"
#include "SLAEBuilder.h"
#include "Boundary.h"
#include "Solvers.h"
#include "Layer.h"

#include <string>
#include <fstream>

using namespace std;

void ThreeInitialVectorsWithMass(
	SLAEBuilder& builder,
	Matrix& A,
	vector<double>& b,
	vector<double>& t,
	map<int, Point>& pointsMap,
	vector<Edge>& bound1,
	vector<Edge>& bound2,
	vector<vector<double>*>& Qs
)
{
	vector<double>* q0 = new vector<double>(A.N);
	vector<double>* q1 = new vector<double>(A.N);
	vector<double>* q2 = new vector<double>(A.N);

	// Set parameters for initial vector q0 =========================================================
	vector<double> lambda = { 0 };
	vector<double> gamma = { 1 };
	function<double(double, double, double)> f = [](double x, double y, double t) { return t * x * x * y; };

	builder.SetLambda(&lambda);
	builder.SetGamma(&gamma);
	builder.SetF(&f);
	// ==============================================================================================

	// Get initial vector q0 ========================================================================
	builder.Build(A, b, t[0]);
	Boundary1(A, b, bound1, pointsMap, t[0]);

	Solvers::LOS(A, *q0, b);
	Qs.push_back(q0);

	A.Clear();
	fill(b.begin(), b.end(), 0.0);
	// ==============================================================================================

	// Get initial vector q1 ========================================================================
	builder.Build(A, b, t[1]);
	Boundary1(A, b, bound1, pointsMap, t[1]);

	Solvers::LOS(A, *q1, b);
	Qs.push_back(q1);

	A.Clear();
	fill(b.begin(), b.end(), 0.0);
	// ==============================================================================================

	// Get initial vector q2 ========================================================================
	builder.Build(A, b, t[2]);
	Boundary1(A, b, bound1, pointsMap, t[2]);

	Solvers::LOS(A, *q2, b);
	Qs.push_back(q2);

	A.Clear();
	fill(b.begin(), b.end(), 0.0);
	// ==============================================================================================
}

void ThreeInitialVectorsWithLayers(
	SLAEBuilder& builder,
	Matrix& A,
	vector<double>& b,
	vector<double>& t,
	map<int, Point>& pointsMap,
	vector<Edge>& bound1,
	vector<Edge>& bound2,
	vector<vector<double>*>& Qs
)
{
	vector<double>* q0 = new vector<double>(A.N);
	vector<double>* q1 = new vector<double>(A.N);
	vector<double>* q2 = new vector<double>(A.N);

	// Set parameters for initial vector q0 =========================================================
	vector<double> lambda = { 0 };
	vector<double> gamma = { 1 };
	function<double(double, double, double)> f = [](double x, double y, double t) { return t * t; };

	builder.SetLambda(&lambda);
	builder.SetGamma(&gamma);
	builder.SetF(&f);
	// ==============================================================================================

	// Get initial vector q0 ========================================================================
	builder.Build(A, b, t[0]);

	Boundary1(A, b, bound1, pointsMap, t[0]);

	Solvers::LOS(A, *q0, b);
	Qs.push_back(q0);

	A.Clear();
	fill(b.begin(), b.end(), 0.0);
	// ==============================================================================================

	// Set parameters for initial vector q0 =========================================================
	lambda = { 1 };
	f = [](double x, double y, double t) { return 2 * t; };

	builder.SetLambda(&lambda);
	builder.SetF(&f);
	// ==============================================================================================

	// Solve two-layer scheme for q1 ================================================================
	Layer* twoLayer = new TwoLayer();

	twoLayer->SetQ({ Qs[0] });
	twoLayer->SetT({ t[0], t[1] });

	builder.SetLayer(twoLayer);
	builder.Build(A, b, t[1]);
	Boundary2(A, b, bound2, pointsMap, lambda[0], t[1]);
	Boundary1(A, b, bound1, pointsMap, t[1]);

	Solvers::LOS(A, *q1, b);
	Qs.push_back(q1);

	fill(b.begin(), b.end(), 0.0);
	A.Clear();
	// ==============================================================================================

	// Solve three-layer scheme for q2 ==============================================================
	Layer* threeLayer = new ThreeLayer();

	threeLayer->SetQ({ Qs[0], Qs[1] });
	threeLayer->SetT({ t[0], t[1], t[2] });

	builder.SetLayer(threeLayer);
	builder.Build(A, b, t[2]);
	Boundary2(A, b, bound2, pointsMap, lambda[0], t[2]);
	Boundary1(A, b, bound1, pointsMap, t[2]);

	Solvers::LOS(A, *q2, b);
	Qs.push_back(q2);

	fill(b.begin(), b.end(), 0.0);
	A.Clear();
	// ==============================================================================================
}

int main()
{
	// Build mesh and time mesh =====================================================================	
	int nodeCount;
	vector<Point> points;
	vector<FiniteElement> elements;

	InputGrid(nodeCount, "points.txt", points, "triangles.txt", elements);
	MeshBuilder mesh(nodeCount);
	mesh.Build(elements);
	nodeCount = mesh.nodeCount;

	vector<Interval> time_intervals;
	InputTime(time_intervals, "time.txt");
	TimeMeshBuilder time(time_intervals);

	vector<double> t = time.getMesh();
	// ==============================================================================================

	// Build boundary conditions ====================================================================
	vector<Edge> bound1, bound2;
	InputBound(bound1, "boundary1.txt");
	mesh.BuildBoundary(bound1);
	InputBound(bound2, "boundary2.txt");
	mesh.BuildBoundary(bound2);
	// ==============================================================================================

	// Build matrix portrait ========================================================================
	Matrix A;
	PortraitBuilder portrait(nodeCount, elements);
	portrait.Build(A);
	// ==============================================================================================

	// SLAE Builder =================================================================================
	SLAEBuilder builder(points, elements);
	vector<double> b(nodeCount);
	vector<double>* q0 = new vector<double>(nodeCount);
	vector<double>* q1 = new vector<double>(nodeCount);
	vector<double>* q2 = new vector<double>(nodeCount);
	vector<vector<double>*> Qs;
	// ==============================================================================================

	// Point map ====================================================================================
	map<int, Point> pointsMap;
	CreatePointMap(elements, points, pointsMap);
	// ==============================================================================================


	ThreeInitialVectorsWithLayers(builder, A, b, t, pointsMap, bound1, bound2, Qs);

	// Set parameters for four-layer scheme =========================================================
	vector<double> lambda = { 1. };
	vector<double> gamma = { 1 };
	function<double(double, double, double)> f = [](double x, double y, double t) { return 2 * t; };
	builder.SetLambda(&lambda);
	builder.SetGamma(&gamma);
	builder.SetF(&f);
	// ==============================================================================================

	// Loop for four-layer scheme ===================================================================
	Layer* fourLayer = new FourLayer();
	builder.SetLayer(fourLayer);
	vector<double>* q;
	for (int i = 3; i < t.size(); i++)
	{
		fourLayer->SetQ({ Qs[i - 3], Qs[i - 2], Qs[i - 1] });
		fourLayer->SetT({ t[i - 3], t[i - 2], t[i - 1], t[i] });

		builder.Build(A, b, t[i]);
		Boundary2(A, b, bound2, pointsMap, lambda[0], t[i]);
		Boundary1(A, b, bound1, pointsMap, t[i]);

		q = new vector<double>(nodeCount);
		Solvers::LOS(A, *q, b);
		Qs.push_back(q);

		fill(b.begin(), b.end(), 0.0);
		A.Clear();
	}
	// ==============================================================================================

	// Output solution ==============================================================================
	for (int i = 1; i < Qs.size(); i++)
	{
		for (int j = 0; j < (*Qs[i]).size(); j++)
			printf("%.10f\n", (*Qs[i])[j]);

		printf("\n");
	}
	// ==============================================================================================

	// Calculation of a function value at an arbitrary point ========================================
	double p1 = 1.0 / 3, p2 = 2.0 / 3;
	vector<double> L = Ls(Point(p1, p2), pointsMap[0], pointsMap[1], pointsMap[2]);
	FiniteElement e = elements[0];

	double res = 0.0;
	vector<double> last = *Qs.back();

	for (int i = 0; i < basisSize; i++)
		res += basisValue(i, L) * last[e.verts[i]];

	printf("Function value at an point (%.5f, %.5f) is %.10f\n\n", p1, p2, res);
	// ==============================================================================================

	system("pause");
	return 0;
}