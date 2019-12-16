#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <set>
#include <functional>

using namespace std;

struct DerComp
{
	int gradNo;
	double coeff;
	int v1, v2, v3;
};

struct PsiComp
{
	double coeff;
	int v1, v2, v3;
};

struct LocalComp
{
	LocalComp(int grad1, int grad2, double coeff) : grad1(grad1), grad2(grad2), coeff(coeff) {}
	int grad1, grad2;
	double coeff;
};

struct Vertex
{
	int vertex, valueNo;
};

struct Point
{
	Point(double x, double y) : x(x), y(y) {}
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

struct Edge
{
	int v1, v2, v3, v4;
	int thetaNo;
};

struct Triangle
{
	vector<int> verts;
	int materialNo;
};

struct Matrix
{
	int N;
	vector<int> IA, JA;
	vector<double> DI, AL;
};

typedef vector<vector<vector<LocalComp>>> GPattern;
typedef vector<vector<double>> MPattern;
typedef vector<vector<double>> LocalMatrix;
typedef vector<set<int>> MatrixPortrait;
typedef vector<vector<int>> EdgeMatrix;

class FEM
{
public:
	FEM(string pointsFile, string trianglesFile, string bounds1File, string bounds2File);

private:
	void InputGrid(string pointsFile, string trianglesFile);
	void InputBound(string bounds1File, string bounds2File);
	void AllocateMemory();
	void VertexNumbering();
	void Portrait();
	void CreateLocalPattern();
	void BuildLocalMatrix(Triangle& t);
	void BuildLocalB(Triangle& t);
	void BuildGlobal();
	void AddLocalToGlobal(Triangle& t);

	int triangleCount;
	int nodeCount;
	int edgeCount;
	int vertexCount;

	vector<Point> points;
	vector<Triangle> triangles;
	vector<Vertex> bound1;
	vector<Edge> bound2;

	EdgeMatrix edgeMatrix;
	MPattern mPattern;
	GPattern gPattern;

	vector<double> alpha;
	vector<int> materials;
	LocalMatrix localMatrix;
	Matrix globalMatrix;
	vector<double> localB;
	vector<double> globalB;

	vector<Point> CalculateCoords(Triangle& t);
	double Alpha(Triangle& t);
	double Det(Triangle& t);
	int Factorial(int N);
};