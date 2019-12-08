#include <iostream>
#include <fstream>
#include <iomanip>
#include <math.h>
#include <vector>
#include <map>
#include <utility>
#include <set>

using namespace std;

int num_triangles, global_number;
double gamma, lambda, L1, L2, L3;
vector<double>  material;
vector<vector<double>> point;
vector<vector<int>>  triangle;
vector<map<int, int>> matr;
vector<set<int>> S;
vector<double> alpha_el;
vector<int> ig, jg;
vector<vector<double>> G, MPattern;
vector<vector<vector<local_el>>> GPattern;
vector<vector<double>> global(global_number);
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
vector<double> elems = { 4.5, 3, 0, 0,          -4.5, 2, 0, 0,       1, 1, 0, 0,
						 4.5, 0, 3, 0,          -4.5, 0, 2, 0,       1, 0, 1, 0,
						 4.5, 0, 0, 3,          -4.5, 0, 0, 2,       1, 0, 0, 1,
						 13.5, 2, 1, 0,         -4.5, 1, 1, 0,       13.5, 1, 2, 0,         -4.5, 1, 1, 0,
						 13.5, 0, 2, 1,         -4.5, 0, 1, 1,       13.5, 0, 1, 2,         -4.5, 0, 1, 1,
						 13.5, 1, 0, 2,         -4.5, 1, 0, 1,       13.5, 2, 0, 1,         -4.5, 1, 0, 1,      27, 1, 1, 1 };

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

vector<vector<der_el>> derivatives;
vector<vector<psi_el>> psi;

void Input()
{
	int num_points;
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
}

void InputDerivatives()
{
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
		for (int j = 0; j < 4; j++)
		{
			psi[i][j].coeff = elems[index++];
			psi[i][j].v1 = elems[index++];
			psi[i][j].v2 = elems[index++];
			psi[i][j].v3 = elems[index++];;
		}
	}

	for (int j = 0; j < 3; j++)
	{
		psi[9][j].coeff = elems[index++];
		psi[9][j].v1 = elems[index++];
		psi[9][j].v2 = elems[index++];
		psi[9][j].v3 = elems[index++];
	}
}

double Alpha(vector<int>& t)
{
	double x1 = point[t[0]][0];
	double y1 = point[t[0]][1];
	double x2 = point[t[1]][0];
	double y2 = point[t[1]][1];
	double x3 = point[t[2]][0];
	double y3 = point[t[2]][1];

	double D = fabs((x2 - x1) * (y3 - y1) - (x3 - x1) * (y2 - y1));

	alpha_el[0] = (y2 - y3) / D;
	alpha_el[1] = (x3 - x2) / D;
	alpha_el[2] = (y3 - y1) / D;
	alpha_el[3] = (x1 - x3) / D;
	alpha_el[4] = (y1 - y2) / D;
	alpha_el[5] = (x2 - x1) / D;

	return D;
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
	matr.resize(global_number);
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
				if (f)swap(a, b);

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

void Portrait()
{
	S.resize(global_number);

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

				S[a].insert(b);
			}
	}

	ig.resize(global_number + 1);
	ig[0] = ig[1] = 0;
	for (int i = 2; i <= global_number; i++)
	{
		int col = ig[i - 1];
		ig[i] = col + S[i].size();
	}
	jg.resize(ig[global_number]);
	for (int i = 2, k = 0; i <= global_number; i++)
	{
		int ig0 = ig[i];
		int ig1 = ig[i + 1];
		for (int j : S[i])
		{
			jg[k] = j;
			k++;
		}
	}
}

void LocalPattern()
{
	GPattern.resize(10);
	for (int i = 0; i < 10; i++)
		GPattern[i].resize(10);

	for (int i = 0; i < 10; i++)
		for (int j = 0; j < 10; j++)
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
			for (int i = 0; i < 10; i++)
				for (int j = i; j < 10; j++)
				{
					double coeff = 0;
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

void BuildGlobal()
{
	for (int i = 0; i < global_number; i++)
		global[i].resize(global_number);

	for (int t = 0; t < num_triangles; t++)
	{
		LocalPattern();

		auto tr = triangle[t];
		double D = Alpha(tr);

		vector<pair<double, double>> grads =
		{
		   make_pair(alpha_el[0], alpha_el[1]),
		   make_pair(alpha_el[2], alpha_el[3]),
		   make_pair(alpha_el[4], alpha_el[5])
		};

		for (int i = 0; i < 10; i++)
		{
			for (int j = 0; j < i + 1; j++)
			{
				for (auto local : GPattern[i][j])
				{
					double scalGrad = grads[local.g1].first * grads[local.g2].first + grads[local.g1].second * grads[local.g2].second;
					G[i][j] += local.coeff * scalGrad * D;
				}
				int i_global = tr[i];
				int j_global = tr[j];

				global[i_global][j_global] += MPattern[i][j] * D + G[i][j];
			}
		}
	}
}

int main()
{
	Input();
	Global_Numbering();
	Portrait();

	return 0;
}