#include <iostream>
#include <fstream>
#include <iomanip>
#include <math.h>
#include <vector>
#include <map>
#include <utility>
#include <set>
#include <functional>

using namespace std;

int num_triangles, global_number;
double L1, L2, L3;
vector<int>  material;
vector<vector<double>> point;
vector<vector<int>> triangle;
vector<double> alpha_el(6);
vector<int> ig, jg;
vector<int> di, al;
vector<vector<double>> G(10), MPattern(10);
vector<vector<double>> global(global_number);
vector<map<int, int>> matr(global_number);
vector<double> B;
vector<double> ders = { 0, 13.5, 2, 0, 0,      0, -9, 1, 0, 0,      0, 1, 0, 0, 0,
						1, 13.5, 0, 2, 0,      1, -9, 0, 2, 0,      1, 1, 0, 0, 0,
						2, 13.5, 0, 0, 3,      2, -9, 0, 0, 3,      2, 1, 0, 0, 0,
						0, 27, 1, 1, 0,        0, -4.5, 0, 1, 0,    1, 13.5, 2, 0, 0,      1, -4.5, 1, 0, 0,
						0, 13.5, 0, 2, 0,      0, -4.5, 0, 1, 0,    1, 27, 1, 1, 0,        1, -4.5, 1, 0, 0,
						1, 27, 0, 1, 1,        1, -4.5, 0, 0, 1,    2, 13.5, 0, 2, 0,      2, -4.5, 0, 1, 0,
						1, 13.5, 0, 0, 2,      1, -4.5, 0, 0, 1,    2, 27, 0, 0, 1,        2, -4.5, 0, 1, 0,
						0, 13.5, 0, 0, 2,      0, -4.5, 0, 0, 1,    2, 27, 1, 0, 1,        2, -4.5, 1, 0, 0,
						0, 27, 1, 0, 1,        0, -4.5, 0, 0, 1,    2, 13.5, 2, 0, 0,      2, -4.5, 1, 0, 0,
						0, 27, 0, 1, 1,        1, 27, 1, 0, 1,      2, 27, 1, 1, 0 };

//vector<int> ders = {    0, 1, 0, 0, 0,         1, 1, 0, 0, 0,       2, 1, 0, 0, 0,
//                        0, 1, 0, 1, 0,         1, 1, 1, 0, 0,
//                        0, 1, 0, 0, 1,         2, 1, 1, 0, 0,
//                        1, 1, 0, 0, 1,         2, 1, 0, 1, 0,
//                        0, 2, 1, 1, 0,         1, 1, 2, 0, 0,       0, -1, 0, 2, 0,        1, -2, 1, 1, 0,
//                        0, 2, 1, 0, 1,         2, 1, 2, 0, 0,       0, -1, 0, 0, 2,        2, -2, 1, 0, 1,
//                        1, 2, 0, 1, 1,         2, 1, 0, 2, 0,       1, -1, 0, 0, 2,        2, -2, 0, 1, 1,
//                        0, 1, 0, 1, 1,         1, 1, 1, 0, 1,       2, 1, 1, 1, 0 };

vector<double> elems = { 4.5, 3, 0, 0,         -4.5, 2, 0, 0,       1, 1, 0, 0,
						 4.5, 0, 3, 0,         -4.5, 0, 2, 0,       1, 0, 1, 0,
						 4.5, 0, 0, 3,         -4.5, 0, 0, 2,       1, 0, 0, 1,
						 13.5, 2, 1, 0,        -4.5, 1, 1, 0,       13.5, 1, 2, 0,         -4.5, 1, 1, 0,
						 13.5, 0, 2, 1,        -4.5, 0, 1, 1,       13.5, 0, 1, 2,         -4.5, 0, 1, 1,
						 13.5, 1, 0, 2,        -4.5, 1, 0, 1,       13.5, 2, 0, 1,         -4.5, 1, 0, 1,      27, 1, 1, 1 };

//vector<int> elems = {    1, 1, 0, 0,         1, 0, 1, 0,          1, 0, 0, 1,
//                         1, 1, 1, 0,         1, 1, 0, 1,          1, 0, 1, 1,    
//                         1, 2, 1, 0,         -1, 1, 2, 0,
//                         1, 2, 0, 1,         -1, 1, 0, 2,
//                         1, 0, 2, 1,         -1, 0, 1, 2,
//                         1, 1, 1, 1 };

vector<double> gamma = { 1, 1, 2 };
vector<double> lambda = { 1, 1, 2 };
vector<double> u_value = { 1, 1, 2 };
vector<double> theta_value = { 1, 1, 2 };
vector<function<double(double, double)>> f = { [](double x, double y) { return x + y; },
											   [](double x, double y) { return x * x; },
											   [](double x, double y) { return y * y; } };


struct der_el
{
	int gradNo;
	double coeff;
	int v1;
	int v2;
	int v3;
};

struct psi_el
{
	double coeff;
	int v1;
	int v2;
	int v3;
};

struct local_el
{
	local_el(int g1, int g2, double coeff) : g1(g1), g2(g2), coeff(coeff) {};
	int g1;
	int g2;
	double coeff;
};

struct point_str
{
	point_str(double x, double y) : x(x), y(y) {};
	double x;
	double y;
};

struct edge
{
	edge(int ver1, int ver2, int ver3, int ver4) : ver1(ver1), ver2(ver2), ver3(ver3), ver4(ver4), theta_num(theta_num) {};
	int ver1;
	int ver2;
	int ver3;
	int ver4;
	int theta_num;
};

struct vertex
{
	vertex(int v, int u_num) : v(v), u_num(u_num) {};
	int v;
	int u_num;
};

vector<vector<der_el>> derivatives;
vector<vector<psi_el>> psi;
vector<edge> bound2;
vector<vertex> bound1;
vector<vector<vector<local_el>>> GPattern(10);

void Input()
{
	int num_points, num_edges, num_vertex;
	ifstream in;
	in.open("points.txt");
	in >> num_points;
	point.resize(num_points);
	for (int i = 0; i < num_points; i++)
	{
		point[i].resize(2);
		in >> point[i][0] >> point[i][1];
	}
	in.close();
	global_number = num_points;

	in.open("triangles.txt");
	in >> num_triangles;
	triangle.resize(num_triangles);
	material.resize(num_triangles);
	for (int i = 0; i < num_triangles; i++)
	{
		triangle[i].resize(10);
		in >> triangle[i][0]  // x1 y1
			>> triangle[i][1]  // x2 y2
			>> triangle[i][2]  // x3 y3
			>> material[i];
	}
	in.close();

	in.open("boundary_conditions_2.txt");
	in >> num_edges;
	bound2.resize(num_edges);
	for (int i = 0; i < num_edges; i++)
	{
		in >> bound2[i].ver1 >> bound2[i].ver4 >> bound2[i].theta_num;
		auto a = bound2[i].ver1;
		auto b = bound2[i].ver4;
		bool f = a < b;
		if (f) swap(a, b);
		bound2[i].ver2 = matr[a][b];
		bound2[i].ver3 = matr[a][b] + 1;
	}
	in.close();

	in.open("boundary_conditions_1.txt");
	in >> num_vertex;
	bound2.resize(num_vertex);
	for (int i = 0; i < num_vertex; i++)
		in >> bound1[i].v >> bound1[i].u_num;
	in.close();
}

void InputDerivatives()
{
	derivatives.resize(10);
	for (int i = 0; i < 3; i++) derivatives[i].resize(3);
	for (int i = 3; i < 9; i++) derivatives[i].resize(4);
	derivatives[9].resize(3);

	int index = 0;

	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			derivatives[i][j].gradNo = ders[index++];
			derivatives[i][j].coeff = ders[index++];
			derivatives[i][j].v1 = ders[index++];
			derivatives[i][j].v2 = ders[index++];
			derivatives[i][j].v3 = ders[index++];
		}
	}

	for (int i = 3; i < 9; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			derivatives[i][j].gradNo = ders[index++];
			derivatives[i][j].coeff = ders[index++];
			derivatives[i][j].v1 = ders[index++];
			derivatives[i][j].v2 = ders[index++];
			derivatives[i][j].v3 = ders[index++];
		}
	}

	for (int j = 0; j < 3; j++)
	{
		derivatives[9][j].gradNo = ders[index++];
		derivatives[9][j].coeff = ders[index++];
		derivatives[9][j].v1 = ders[index++];
		derivatives[9][j].v2 = ders[index++];
		derivatives[9][j].v3 = ders[index++];
	}
}

void InputFunctions()
{
	psi.resize(10);
	for (int i = 0; i < 3; i++) psi[i].resize(3);
	for (int i = 3; i < 9; i++) psi[i].resize(4);
	psi[9].resize(3);

	int index = 0;

	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			psi[i][j].coeff = elems[index++];
			psi[i][j].v1 = elems[index++];
			psi[i][j].v2 = elems[index++];
			psi[i][j].v3 = elems[index++];
		}
	}

	for (int i = 3; i < 9; i++)
	{
		for (int j = 0; j < 2; j++)
		{
			psi[i][j].coeff = elems[index++];
			psi[i][j].v1 = elems[index++];
			psi[i][j].v2 = elems[index++];
			psi[i][j].v3 = elems[index++];;
		}
	}

	psi[9][0].coeff = elems[index++];
	psi[9][0].v1 = elems[index++];
	psi[9][0].v2 = elems[index++];
	psi[9][0].v3 = elems[index++];
}

double Alpha(vector<int>& t)
{
	double x1 = point[t[0]][0];
	double y1 = point[t[0]][1];
	double x2 = point[t[1]][0];
	double y2 = point[t[1]][1];
	double x3 = point[t[2]][0];
	double y3 = point[t[2]][1];

	double D = det(t);

	alpha_el[0] = (y2 - y3) / D;
	alpha_el[1] = (x3 - x2) / D;
	alpha_el[2] = (y3 - y1) / D;
	alpha_el[3] = (x1 - x3) / D;
	alpha_el[4] = (y1 - y2) / D;
	alpha_el[5] = (x2 - x1) / D;
}

double det(vector<int>& t)
{
	double x1 = point[t[0]][0];
	double y1 = point[t[0]][1];
	double x2 = point[t[1]][0];
	double y2 = point[t[1]][1];
	double x3 = point[t[2]][0];
	double y3 = point[t[2]][1];

	double D = (x2 - x1) * (y3 - y1) - (x3 - x1) * (y2 - y1);

	return fabs(D);
}

//void L_coord(int i, double x, double y)
//{
//   double x1 = point[triangle[i][0]][0], x2 = point[triangle[i][1]][0], x3 = point[triangle[i][2]][0];
//   double y1 = point[triangle[i][0]][1], y2 = point[triangle[i][1]][1], y3 = point[triangle[i][2]][1];
//   double D = (x2 - x1) * (y3 - y1) - (x3 - x1) * (y2 - y1);
//   L1 = fabs(((x3 - x2) * (y - y2) - (x - x2) * (y3 - y2)) / D);
//   L2 = fabs(((x1 - x3) * (y - y3) - (x - x3) * (y1 - y3)) / D);
//   L3 = fabs(((x2 - x1) * (y - y1) - (x - x1) * (y2 - y1)) / D);
//}

int Fact(int N)
{
	int res = 1;
	for (int i = 2; i <= N; i++) res *= i;
	return res;
}

void Global_Numbering()
{
	for (int t = 0; t < num_triangles; t++)
	{
		int index = 3;
		for (int i = 1; i < 3; i++)
		{
			for (int j = 0; j < i; j++, index += 2)
			{
				auto a = triangle[t][i];
				auto b = triangle[t][j];
				bool f = a < b;
				if (f) swap(a, b);

				if (matr[a][b] == 0)
				{
					triangle[t][index] = global_number + (f ? 0 : 1);
					triangle[t][index + 1] = global_number + (f ? 1 : 0);
					matr[a][b] = global_number;
					global_number += 2;
				}
				else
				{
					auto a1 = matr[a][b];
					triangle[t][index] = a1 + (f ? 0 : 1);
					triangle[t][index + 1] = a1 + (f ? 1 : 0);
				}
			}
		}
		triangle[t][index] = global_number;
		global_number++;
	}
}

void Boundary2()
{

}

void Portrait()
{
	vector<set<int>> Connection(global_number);

	for (int t = 0; t < num_triangles; t++)
	{
		int index = 3;
		for (int i = 1; i < 10; i++)
			for (int j = 0; j < i; j++, index += 2)
			{
				auto a = triangle[t][i];
				auto b = triangle[t][j];
				bool f = a < b;
				if (f)swap(a, b);

				Connection[a].insert(b);
			}
	}

	ig.resize(global_number + 1);
	ig[0] = ig[1] = 0;
	for (int i = 2; i <= global_number; i++)
	{
		int col = ig[i - 1];
		ig[i] = col + Connection[i - 1].size();
	}
	jg.resize(ig[global_number]);
	for (int i = 1, k = 0; i < global_number; i++)
	{
		for (int j : Connection[i])
		{
			jg[k] = j;
			k++;
		}
	}
	di.resize(global_number);
	al.resize(ig[global_number]);
}

void Memory()
{
	B.resize(global_number);
	for (int i = 0; i < 10; i++)
	{
		G[i].resize(i + 1);
		GPattern[i].resize(i + 1);
		MPattern[i].resize(i + 1);
	}
}

void LocalPattern()
{
	for (int i = 0; i < 10; i++)
	{
		for (int j = 0; j <= i; j++)
		{
			for (auto x : derivatives[i])
				for (auto y : derivatives[j])
				{
					int v1 = x.v1 + y.v1;
					int v2 = x.v2 + y.v2;
					int v3 = x.v3 + y.v3;


					double coeff = x.coeff * y.coeff * Fact(v1) * Fact(v2) * Fact(v3) / Fact(v1 + v2 + v3 + 2);
					GPattern[i][j].push_back(local_el(x.gradNo, y.gradNo, coeff));
				}

			double sum = 0;
			for (auto x : psi[i])
				for (auto y : psi[j])
				{
					int v1 = x.v1 + y.v1;
					int v2 = x.v2 + y.v2;
					int v3 = x.v3 + y.v3;

					sum += x.coeff * y.coeff * Fact(v1) * Fact(v2) * Fact(v3) / Fact(v1 + v2 + v3 + 2);
				}
			MPattern[i][j] = sum;
		}
	}
}

void BuildLocal(int t, vector<vector<double>>& Local)
{
	auto tr = triangle[t];
	double D = det(tr);

	vector<pair<double, double>> grads =
	{
	   make_pair(alpha_el[0], alpha_el[1]),
	   make_pair(alpha_el[2], alpha_el[3]),
	   make_pair(alpha_el[4], alpha_el[5])
	};

	for (int i = 0; i < 10; i++)
	{
		G[i].resize(i + 1);
		for (int j = 0; j <= i; j++)
		{
			for (auto local : GPattern[i][j])
			{
				double scalGrad = grads[local.g1].first * grads[local.g2].first + grads[local.g1].second * grads[local.g2].second;
				G[i][j] += local.coeff * scalGrad;
			}
			Local[i][j] = (gamma[material[t]] * MPattern[i][j] + lambda[material[t]] * G[i][j]) * D;
		}
	}
}

void AddToGlobal(int t, vector<vector<double>>& Local, vector<double>& b)
{
	auto tr = triangle[t];
	for (int i = 0; i < 10; i++)
	{
		di[tr[i]] += Local[i][i];
		B[tr[i]] += b[i];

		int beg = ig[tr[i]];
		for (int j = 0; j < i - 1; j++, beg++)
		{
			int end = ig[tr[i] + 1] - 1;
			while (jg[beg] != tr[j])
			{
				int ind = (beg + end) / 2 + 1;
				if (jg[ind] <= tr[j]) beg = ind;
				else end = ind;
			}
			al[beg] += Local[i][j];
		}
	}
}

void BuildGlobal()
{
	vector<vector<double>> Local(10);
	vector<double> b(10);
	for (int i = 0; i < 10; i++) Local[i].resize(i + 1);
	LocalPattern();
	for (int t = 0; t < num_triangles; t++)
	{
		auto tr = triangle[t];
		BuildLocal(t, Local);
		b = BuildLocalB(tr, material[t]);
		AddToGlobal(t, Local, b);
		Local.clear();
		b.clear();
	}
}

vector<double> BuildLocalB(vector<int>& t, int m)
{
	vector<point_str> coords;
	double D = det(t);
	coords.push_back(point_str(point[t[0]][0], point[t[0]][1]));
	coords.push_back(point_str(point[t[1]][0], point[t[1]][1]));
	coords.push_back(point_str(point[t[2]][0], point[t[2]][1]));
	coords.push_back(point_str((coords[0].x + coords[1].x) / 3., (coords[0].y + coords[1].y) / 3.));
	coords.push_back(point_str(2 * (coords[0].x + coords[1].x) / 3., 2 * (coords[0].y + coords[1].y) / 3.));
	coords.push_back(point_str((coords[1].x + coords[2].x) / 3., (coords[1].y + coords[2].y) / 3.));
	coords.push_back(point_str(2 * (coords[1].x + coords[2].x) / 3., 2 * (coords[1].y + coords[2].y) / 3.));
	coords.push_back(point_str((coords[0].x + coords[2].x) / 3., (coords[0].y + coords[2].y) / 3.));
	coords.push_back(point_str(2 * (coords[0].x + coords[2].x) / 3., 2 * (coords[0].y + coords[2].y) / 3.));
	coords.push_back(point_str((coords[0].x + coords[1].x + coords[2].x) / 3., (coords[0].y + coords[1].y + coords[2].y) / 3.));

	vector<double> b(coords.size(), 0.0);
	vector<double> temp(coords.size(), 0.0);

	for (int i = 0; i < 10; i++)
		temp[i] = f[m](coords[i].x, coords[i].y);

	for (int i = 0; i < 10; i++)
		for (int j = 0; j < 10; j++)
			b[i] += temp[j] * MPattern[i][j] * D;

	return b;
}

int main()
{
	Input();
	InputDerivatives();
	InputFunctions();
	Global_Numbering();
	Portrait();
	Memory();
	BuildGlobal();

	return 0;
}