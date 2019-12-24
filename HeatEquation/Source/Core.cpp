#include "Core.h"

#include <iostream>
#include <iomanip>



vector<vector<DerComp>> ders = {
	{ { 0, 13.5, 2, 0, 0 },			{ 0, -9, 1, 0, 0 },			{ 0, 1, 0, 0, 0 } },
	{ { 1, 13.5, 0, 2, 0 },			{ 1, -9, 0, 1, 0 },			{ 1, 1, 0, 0, 0 } },
	{ { 2, 13.5, 0, 0, 2 },			{ 2, -9, 0, 0, 1 },			{ 2, 1, 0, 0, 0 } },
	{ { 0, 27, 1, 1, 0 },			{ 0, -4.5, 0, 1, 0 },		{ 1, 13.5, 2, 0, 0 },		{ 1, -4.5, 1, 0, 0 } },
	{ { 0, 13.5, 0, 2, 0 },			{ 0, -4.5, 0, 1, 0 },		{ 1, 27, 1, 1, 0 },			{ 1, -4.5, 1, 0, 0 } },
	{ { 1, 27, 0, 1, 1 },			{ 1, -4.5, 0, 0, 1 },		{ 2, 13.5, 0, 2, 0 },		{ 2, -4.5, 0, 1, 0 } },
	{ { 1, 13.5, 0, 0, 2 },			{ 1, -4.5, 0, 0, 1 },		{ 2, 27, 0, 1, 1 },			{ 2, -4.5, 0, 1, 0 } },
	{ { 0, 13.5, 0, 0, 2 },			{ 0, -4.5, 0, 0, 1 },		{ 2, 27, 1, 0, 1 },			{ 2, -4.5, 1, 0, 0 } },
	{ { 0, 27, 1, 0, 1 },			{ 0, -4.5, 0, 0, 1 },		{ 2, 13.5, 2, 0, 0 },		{ 2, -4.5, 1, 0, 0 } },
	{ { 0, 27, 0, 1, 1 },			{ 1, 27, 1, 0, 1 },			{ 2, 27, 1, 1, 0 } }
};

vector<vector<PsiComp>> basis = {
	{ { 4.5, 3, 0, 0 },         { -4.5, 2, 0, 0 },       { 1, 1, 0, 0 } },
	{ { 4.5, 0, 3, 0 },         { -4.5, 0, 2, 0 },       { 1, 0, 1, 0 } },
	{ { 4.5, 0, 0, 3 },         { -4.5, 0, 0, 2 },       { 1, 0, 0, 1 } },
	{ { 13.5, 2, 1, 0 },        { -4.5, 1, 1, 0 } },
	{ { 13.5, 1, 2, 0 },        { -4.5, 1, 1, 0 } },
	{ { 13.5, 0, 2, 1 },        { -4.5, 0, 1, 1 } },
	{ { 13.5, 0, 1, 2 },		{ -4.5, 0, 1, 1 } },
	{ { 13.5, 1, 0, 2 },        { -4.5, 1, 0, 1 } },
	{ { 13.5, 2, 0, 1 },        { -4.5, 1, 0, 1 } },
	{ { 27, 1, 1, 1 } }
};

vector<function<double(double, double)>> f = {
	[](double x, double y) { return 0.0; },
	[](double x, double y) { return x * x; },
	[](double x, double y) { return y * y; }
};

vector<double> gamma = { 0, 0, 0 };
vector<double> lambda = { 1, 1, 1 };
vector<double> uValue = { -1, 1 };
vector<double> thetaValue = { 0, 0, 0 };
vector<double> edgeBasisValues = { 0.125, 0.125, 0.375, 0.375 };

FEM::FEM(string pointsFile, string trianglesFile, string bounds1File, string bounds2File)
{
	InputGrid(pointsFile, trianglesFile);
	AllocateMemory();
	VertexNumbering();
	Portrait();
	InputBound(bounds1File, bounds2File);
	CreateLocalPattern();
	BuildGlobal();
	//Boundary2();
	Boundary1();
	Compute();
	Output();
}

void FEM::InputGrid(string pointsFile, string trianglesFile)
{
	ifstream in(pointsFile);
	in >> nodeCount;
	points.resize(nodeCount);
	for (int i = 0; i < nodeCount; i++)
		in >> points[i].x >> points[i].y;
	in.close();

	in.open(trianglesFile);
	in >> triangleCount;
	triangles.resize(triangleCount);
	for (int i = 0; i < triangleCount; i++)
	{
		triangles[i].verts.resize(10);
		in >> triangles[i].verts[0]
			>> triangles[i].verts[1]
			>> triangles[i].verts[2]
			>> triangles[i].materialNo;
	}
	in.close();
}

void FEM::InputBound(string bounds1File, string bounds2File)
{
	ifstream in(bounds1File);
	in >> edgeCount1;
	bound1.resize(edgeCount1);
	for (int i = 0; i < edgeCount1; i++)
	{
		in >> bound1[i].v1 >> bound1[i].v4 >> bound1[i].valueNo;
		auto a = bound1[i].v1;
		auto b = bound1[i].v4;

		if (b < a) swap(a, b);

		bound1[i].v2 = edgeMatrix[a][b];
		bound1[i].v3 = edgeMatrix[a][b] + 1;
	}
	in.close();


	in.open(bounds2File);
	in >> edgeCount2;
	bound2.resize(edgeCount2);
	for (int i = 0; i < edgeCount2; i++)
	{
		in >> bound2[i].v1 >> bound2[i].v4 >> bound2[i].valueNo;
		auto a = bound2[i].v1;
		auto b = bound2[i].v4;

		if (a < b) swap(a, b);

		bound2[i].v2 = edgeMatrix[a][b];
		bound2[i].v3 = edgeMatrix[a][b] + 1;
	}
	in.close();
}

void FEM::AllocateMemory()
{
	edgeMatrix.resize(nodeCount);
	for (int i = 0; i < nodeCount; i++)
		edgeMatrix[i].resize(nodeCount);

	localMatrix.resize(10);
	gPattern.resize(10);
	mPattern.resize(10);

	for (int i = 0; i < 10; i++)
	{
		localMatrix[i].resize(10);
		gPattern[i].resize(10);
		mPattern[i].resize(10);
	}

	localB.resize(10);
	alpha.resize(6);
}

void FEM::VertexNumbering()
{
	for (auto& t : triangles)
	{
		int index = 3;
		for (int i = 0; i < 3; i++)
			for (int j = i + 1; j < 3; j++, index += 2)
			{
				int a = t.verts[i];
				int b = t.verts[j];
				bool f = a > b;
				if (f) swap(a, b);

				if (edgeMatrix[a][b] == 0)
				{
					t.verts[index] = nodeCount + (f ? 0 : 1);
					t.verts[index + 1] = nodeCount + (f ? 1 : 0);
					edgeMatrix[a][b] = nodeCount;
					nodeCount += 2;
				}
				else
				{
					int a1 = edgeMatrix[a][b];
					t.verts[index] = a1 + (f ? 0 : 1);
					t.verts[index + 1] = a1 + (f ? 1 : 0);
				}
			}
		t.verts[index] = nodeCount;
		nodeCount++;
	}
}

void FEM::Portrait()
{
	MatrixPortrait connection(nodeCount);

	for (auto& t : triangles)
	{
		for (int i = 1; i < 10; i++)
			for (int j = 0; j < i; j++)
			{
				auto a = t.verts[i];
				auto b = t.verts[j];
				if (a < b) swap(a, b);

				connection[a].insert(b);
			}
	}

	globalMatrix.N = nodeCount;
	globalMatrix.IA.resize(nodeCount + 1);
	int* IA = &globalMatrix.IA[0];
	IA[0] = IA[1] = 0;

	for (int i = 2; i <= nodeCount; i++)
	{
		int col = IA[i - 1];
		IA[i] = col + connection[i - 1].size();
	}

	globalMatrix.JA.resize(IA[nodeCount]);
	int* JA = &globalMatrix.JA[0];
	for (int i = 1, k = 0; i < nodeCount; i++)
		for (int j : connection[i])
		{
			JA[k] = j;
			k++;
		}

	globalMatrix.DI.resize(nodeCount);
	globalMatrix.AL.resize(IA[nodeCount]);
	globalB.resize(nodeCount);
	q.resize(nodeCount);
}

void FEM::CreateLocalPattern()
{
	for (int i = 0; i < 10; i++)
	{
		for (int j = 0; j < 10; j++)
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
			mPattern[i][j] = sum;
		}
	}
}

void FEM::BuildLocalMatrix(Triangle& t)
{
	double D = AbsDet(t);
	Alpha(t);
	vector<Grad> grads =
	{
		{ alpha[0], alpha[1] },
		{ alpha[2], alpha[3] },
		{ alpha[4], alpha[5] }
	};

	for (int i = 0; i < 10; i++)
	{
		for (int j = 0; j < 10; j++)
		{
			double sum = 0;
			for (auto comp : gPattern[i][j])
			{
				double scalGrad = grads[comp.grad1].a1 * grads[comp.grad2].a1 + grads[comp.grad1].a2 * grads[comp.grad2].a2;
				sum += comp.coeff * scalGrad;
			}
			localMatrix[i][j] = (gamma[t.materialNo] * mPattern[i][j] + lambda[t.materialNo] * sum) * D;
		}
	}
}

void FEM::BuildLocalB(Triangle& t)
{
	vector<Point> coords = CalculateCoords(t);
	double D = AbsDet(t);

	fill(localB.begin(), localB.end(), 0.0);
	vector<double> temp(10);
	for (int i = 0; i < 10; i++)
		temp[i] = f[t.materialNo](coords[i].x, coords[i].y);

	for (int i = 0; i < 10; i++)
	{
		for (int j = 0; j < 10; j++)
			localB[i] += temp[j] * mPattern[i][j];

		localB[i] *= D;
	}
}

void FEM::AddLocalToGlobal(Triangle& t)
{
	for (int i = 0; i < 10; i++)
	{
		globalMatrix.DI[t.verts[i]] += localMatrix[i][i];
		globalB[t.verts[i]] += localB[i];
		for (int j = 0; j < i; j++)
		{
			auto a = t.verts[i];
			auto b = t.verts[j];
			if (a < b) swap(a, b);

			auto begin = globalMatrix.JA.begin() + globalMatrix.IA[a];
			if (globalMatrix.IA[a + 1] > globalMatrix.IA[a])
			{
				auto end = globalMatrix.JA.begin() + globalMatrix.IA[a + 1] - 1;
				auto iter = lower_bound(begin, end, b);
				auto index = iter - globalMatrix.JA.begin();
				globalMatrix.AL[index] += localMatrix[i][j];
			}
		}
	}

	/*-----------------------------------Tight format---------------------------------------*/
	for (int i = 0; i < 10; i++)
	{
		//globalB[t.verts[i]] += localB[i];
		for (int j = 0; j < 10; j++)
		{
			A[t.verts[i]][t.verts[j]] += localMatrix[i][j];
		}
	}
}

void FEM::BuildGlobal()
{
	/*---------------------------------For debug-----------------------------------------*/
	A.resize(nodeCount);
	for (int i = 0; i < nodeCount; i++)
	{
		A[i].resize(nodeCount);
	}
	/*---------------------------------For debug-----------------------------------------*/

	for (auto& t : triangles)
	{
		BuildLocalMatrix(t);

		/*---------------------------------For debug-----------------------------------------*/
		vector<double> aa = { -1, -1, 1, -1, -1, -1.0 / 3, 1.0 / 3, -1.0 / 3, 1.0 / 3, -1.0 / 3 };
		vector<double> b(10, 0.0);

		for (int i = 0; i < 10; i++)
		{
			for (int j = 0; j < 10; j++)
			{
				b[i] += aa[j] * localMatrix[i][j];
			}
		}
		/*---------------------------------For debug-----------------------------------------*/

		BuildLocalB(t);
		AddLocalToGlobal(t);
	}

	/*---------------------------------For debug-----------------------------------------*/
	vector<vector<double>> A1;
	A1.resize(nodeCount);
	for (int i = 0; i < nodeCount; i++)
	{
		A1[i].resize(nodeCount);
	}

	for (int i = 0; i < nodeCount; i++)
		A1[i][i] = globalMatrix.DI[i];

	for (int i = 0; i < nodeCount; i++)
	{
		int i0 = globalMatrix.IA[i], i1 = globalMatrix.IA[i + 1];

		for (int k = i0; k < i1; k++)
		{
			int j = globalMatrix.JA[k];
			A1[i][j] = globalMatrix.AL[k];
			A1[j][i] = globalMatrix.AL[k];
		}
	}

	for (int i = 0; i < nodeCount; i++)
	{
		for (int j = 0; j < nodeCount; j++)
			cout << A1[i][j] << "\t";

		cout << endl;
	}
	cout << endl << endl;
	for (int i = 0; i < nodeCount; i++)
	{
		double sum = 0;
		for (int j = 0; j < nodeCount; j++)
			sum += A1[i][j];

		cout << sum << endl;
	}
	cout << endl << endl;
	cout << endl << endl << endl << endl;


	for (int i = 0; i < nodeCount; i++)
	{
		for (int j = 0; j < nodeCount; j++)
			cout << A[i][j] << "\t";

		cout << endl;
	}
	cout << endl << endl;

	for (int i = 0; i < nodeCount; i++)
	{
		double sum = 0;
		for (int j = 0; j < nodeCount; j++)
			sum += A[i][j];

		cout << sum << endl;
	}
	cout << endl << endl;
	/*---------------------------------For debug-----------------------------------------*/
}

void FEM::Boundary1()
{
	for (auto edge : bound1)
	{
		globalMatrix.DI[edge.v1] = 1.0E+50;
		globalB[edge.v1] = 1.0E+50 * uValue[edge.valueNo];
		globalMatrix.DI[edge.v2] = 1.0E+50;
		globalB[edge.v2] = 1.0E+50 * uValue[edge.valueNo];
		globalMatrix.DI[edge.v3] = 1.0E+50;
		globalB[edge.v3] = 1.0E+50 * uValue[edge.valueNo];
		globalMatrix.DI[edge.v4] = 1.0E+50;
		globalB[edge.v4] = 1.0E+50 * uValue[edge.valueNo];
	}
}

void FEM::Boundary2()
{
	vector<double> conditions;
	conditions.resize(4);
	edgeBasisValues.resize(4);

	for (int edge = 0; edge < edgeCount2; edge++)
	{
		fill(conditions.begin(), conditions.end(), 0);
		for (int i = 0; i < 4; i++)
		{
			double h = Distance(points[bound2[edge].v1], points[bound2[edge].v4]);
			conditions[i] = h * edgeBasisValues[i] * thetaValue[bound2[edge].valueNo];
		}

		globalB[bound2[edge].v1] += conditions[0];
		globalB[bound2[edge].v4] += conditions[1];
		globalB[bound2[edge].v2] += conditions[2];
		globalB[bound2[edge].v3] += conditions[3];
	}
}

double FEM::Distance(Point a, Point b)
{
	return sqrt((a.x - b.x) * (a.x - b.x) + (a.y - b.y) * (a.y - b.y));
}

void FEM::Compute()
{
	vector<double> r, z, Az;
	int maxiter = 100000;
	double eps = 1.E-15;
	r.resize(nodeCount);
	z.resize(nodeCount);
	Az.resize(nodeCount);

	/*---------------------------------For debug-----------------------------------------*/
	vector<double> aa = { -1, -1, 1, 1, -1, -1, 1.0 / 3, -1.0 / 3, -1.0 / 3, 1.0 / 3, -1.0 / 3, -1.0 / 3, 1.0 / 3, 1.0, 1, 1.0 / 3 };

	//Multiply(aa, q);
	/*---------------------------------For debug-----------------------------------------*/


	// r = A*x
	Multiply(q, r);

	// r = f - A*x 
	for (int i = 0; i < nodeCount; i++)
	{
		r[i] = globalB[i] - r[i];
		z[i] = r[i];
	}

	double normB = 1;//sqrt(fabs(Scal(globalB, globalB)));
	double scal_rr = Scal(r, r);
	double diff = sqrt(fabs(scal_rr)) / normB;
	int iter = 0;
	for (; iter < maxiter && diff >= eps; iter++)
	{
		Multiply(z, Az);     // Az = A*z
		double a = scal_rr / Scal(Az, z);
		for (int i = 0; i < nodeCount; i++)
		{
			q[i] += a * z[i];
			r[i] -= a * Az[i];
		}

		double b = 1. / scal_rr;
		scal_rr = Scal(r, r);
		b *= scal_rr;
		for (int i = 0; i < nodeCount; i++)
			z[i] = r[i] + b * z[i];

		diff = sqrt(fabs(Scal(r, r))) / normB;
	}
}

vector<Point> FEM::CalculateCoords(Triangle& t)
{
	Point t1 = points[t.verts[0]];
	Point t2 = points[t.verts[1]];
	Point t3 = points[t.verts[2]];

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

void FEM::Alpha(Triangle& t)
{
	double x1 = points[t.verts[0]].x;
	double y1 = points[t.verts[0]].y;
	double x2 = points[t.verts[1]].x;
	double y2 = points[t.verts[1]].y;
	double x3 = points[t.verts[2]].x;
	double y3 = points[t.verts[2]].y;

	double D = Det(t);

	alpha[0] = (y2 - y3) / D;
	alpha[1] = (x3 - x2) / D;
	alpha[2] = (y3 - y1) / D;
	alpha[3] = (x1 - x3) / D;
	alpha[4] = (y1 - y2) / D;
	alpha[5] = (x2 - x1) / D;
}

double FEM::AbsDet(Triangle& t)
{
	double x1 = points[t.verts[0]].x;
	double y1 = points[t.verts[0]].y;
	double x2 = points[t.verts[1]].x;
	double y2 = points[t.verts[1]].y;
	double x3 = points[t.verts[2]].x;
	double y3 = points[t.verts[2]].y;

	return abs((x2 - x1) * (y3 - y1) - (x3 - x1) * (y2 - y1));
}

double FEM::Det(Triangle& t)
{
	double x1 = points[t.verts[0]].x;
	double y1 = points[t.verts[0]].y;
	double x2 = points[t.verts[1]].x;
	double y2 = points[t.verts[1]].y;
	double x3 = points[t.verts[2]].x;
	double y3 = points[t.verts[2]].y;

	return (x2 - x1) * (y3 - y1) - (x3 - x1) * (y2 - y1);
}

int FEM::Factorial(int n)
{
	int res = 1;
	for (int i = 2; i <= n; i++)
		res *= i;
	return res;
}

void FEM::Forward(vector<double>& y, vector<double>& b)
{
	for (int i = 0; i < nodeCount; i++)
	{
		double sum = 0;
		int i0 = globalMatrix.IA[i], i1 = globalMatrix.IA[i + 1];
		for (int k = i0; k < i1; k++)
		{
			int j = globalMatrix.JA[k];
			sum += LLT.AL[k] * y[j];
		}
		y[i] = (b[i] - sum) / LLT.DI[i];
	}
}

void FEM::Backward(vector<double>& x, vector<double>& y)
{
	for (int i = 0; i < nodeCount; i++)
		x[i] = y[i];
	for (int i = nodeCount - 1; i >= 0; i--)
	{
		int i0 = globalMatrix.IA[i], i1 = globalMatrix.IA[i + 1];
		x[i] /= LLT.DI[i];
		for (int k = i0; k < i1; k++)
		{
			int j = globalMatrix.JA[k];
			x[j] -= LLT.AL[k] * x[i];
		}
	}
}

void FEM::Multiply(vector<double>& x, vector<double>& res)
{
	for (int i = 0; i < nodeCount; i++)
	{
		res[i] = globalMatrix.DI[i] * x[i];
		int i0 = globalMatrix.IA[i];
		int i1 = globalMatrix.IA[i + 1];
		for (int k = i0; k < i1; k++)
		{
			int j = globalMatrix.JA[k];
			res[i] += globalMatrix.AL[k] * x[j];
			res[j] += globalMatrix.AL[k] * x[i];
		}
	}
}

double FEM::Scal(vector<double>& x, vector<double>& y)
{
	double scal = 0;
	for (int i = 0; i < nodeCount; i++)
		//if (x[i] < 10e+40 && y[i] < 10e+40)
		scal += x[i] * y[i];

	return scal;
}

void FEM::Output()
{
	for (int i = 0; i < nodeCount; i++)
		printf("%.10f\n", q[i]);
}