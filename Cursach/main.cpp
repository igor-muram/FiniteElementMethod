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
vector<double>  material;
vector<vector<double>> point;
vector<vector<int>>  triangle;
vector<map<int, int>> matr;
vector<set<int>> S;
vector<double> alpha_el;
vector<int> ig, jg;

struct der_el
{
	der_el(int gradNo, double coeff, int v1, int v2, int v3) :gradNo(gradNo), coeff(coeff), v1(v1), v2(v2), v3(v3) {}
	int gradNo;
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
vector<vector<der_el>> phi_struct;

vector<vector<vector<local_el>>> localPattern;
double L1, L2, L3;

void input()
{

	int num_points;
	ifstream points;
	points.open("points.txt");
	points >> num_points;
	point.resize(num_points);
	for (int i = 0; i < num_points; i++)
	{
		point[i].resize(2);
		points >> point[i][0] >> point[i][1];  // x y
	}
	points.close();
	global_number = num_points;


	ifstream triangles;
	triangles.open("triangles.txt");
	triangles >> num_triangles;
	triangle.resize(num_triangles);
	material.resize(num_triangles);
	for (int i = 0; i < num_triangles; i++)
	{
		triangle[i].resize(10);
		triangles >> triangle[i][0]  // x1 y1
			>> triangle[i][1]  // x2 y2
			>> triangle[i][2]; // x3 y3
		triangles >> material[i];
	}
	triangles.close();
}

void alpha(vector<int>& t)
{
	double x1 = point[t[0]][0];
	double y1 = point[t[0]][1];
	double x2 = point[t[1]][0];
	double y2 = point[t[1]][1];
	double x3 = point[t[2]][0];
	double y3 = point[t[2]][1];

	double D = abs((x2 - x1) * (y3 - y1) - (x3 - x1) * (y2 - y1));
	alpha_el[0] = (y2 - y3) / D;
	alpha_el[1] = (x3 - x2) / D;
	alpha_el[2] = (y3 - y1) / D;
	alpha_el[3] = (x1 - x3) / D;
	alpha_el[4] = (y1 - y2) / D;
	alpha_el[5] = (x2 - x1) / D;

}

void derivatives_a(int i, int j, int k, bool truex)///////кому не нравится перепишете 
{
	int index = 1;
	if (truex)index--;
	alpha(i, j, k);
	for (int i = 0; i < 3; i++)
	{
		derivatives[i][0].a = 27. / 2 * alpha_el[2 * i];
		derivatives[i][1].a = -9 * alpha_el[2 * i];
		derivatives[i][2].a = alpha_el[2 * i];
	}
	derivatives[3][0].a = 27 * alpha_el[0 + index]; derivatives[3][1].a = 27. / 2 * alpha_el[2 + index];
	derivatives[4][0].a = 27 * alpha_el[2 + index]; derivatives[4][1].a = 27. / 2 * alpha_el[0 + index];
	derivatives[5][0].a = 27 * alpha_el[2 + index]; derivatives[5][1].a = 27. / 2 * alpha_el[4 + index];
	derivatives[6][0].a = 27 * alpha_el[4 + index]; derivatives[6][1].a = 27. / 2 * alpha_el[2 + index];
	derivatives[7][0].a = 27 * alpha_el[4 + index]; derivatives[7][1].a = 27. / 2 * alpha_el[0 + index];
	derivatives[8][0].a = 27 * alpha_el[0 + index]; derivatives[8][1].a = 27. / 2 * alpha_el[4 + index];

	derivatives[3][2].a = -9. / 2 * alpha_el[0 + index]; derivatives[3][3].a = -9. / 2. * alpha_el[2 + index];
	derivatives[4][2].a = -9. / 2 * alpha_el[2 + index]; derivatives[4][3].a = -9. / 2. * alpha_el[0 + index];
	derivatives[5][2].a = -9. / 2 * alpha_el[2 + index]; derivatives[5][3].a = -9. / 2. * alpha_el[4 + index];
	derivatives[6][2].a = -9. / 2 * alpha_el[4 + index]; derivatives[6][3].a = -9. / 2. * alpha_el[2 + index];
	derivatives[7][2].a = -9. / 2 * alpha_el[4 + index]; derivatives[7][3].a = -9. / 2. * alpha_el[0 + index];
	derivatives[8][2].a = -9. / 2 * alpha_el[0 + index]; derivatives[8][3].a = -9. / 2. * alpha_el[4 + index];


	derivatives[9][0].a = 27 * alpha_el[0]; derivatives[9][1].a = 27 * alpha_el[2]; derivatives[9][2].a = 27 * alpha_el[4];
}

void L_coord(int i, double x, double y)
{
	double x1 = point[triangle[i][0]][0], x2 = point[triangle[i][1]][0], x3 = point[triangle[i][2]][0];
	double y1 = point[triangle[i][0]][1], y2 = point[triangle[i][1]][1], y3 = point[triangle[i][2]][1];
	double D = (x2 - x1) * (y3 - y1) - (x3 - x1) * (y2 - y1);
	L1 = fabs(((x3 - x2) * (y - y2) - (x - x2) * (y3 - y2)) / D);
	L2 = fabs(((x1 - x3) * (y - y3) - (x - x3) * (y1 - y3)) / D);
	L3 = fabs(((x2 - x1) * (y - y1) - (x - x1) * (y2 - y1)) / D);
}

double detD(vector<int>& t)
{
	double x1 = point[t[0]][0];
	double y1 = point[t[0]][1];
	double x2 = point[t[1]][0];
	double y2 = point[t[1]][1];
	double x3 = point[t[2]][0];
	double y3 = point[t[2]][1];

	double D = (x2 - x1) * (y3 - y1) - (x3 - x1) * (y2 - y1);

	return D;
}

void basics_function()
{
	ifstream phi;
	phi.open("phi.txt");
	for (int i = 0; i < 3; i++) phi_struct[i].resize(3);
	for (int i = 3; i < 9; i++) phi_struct[i].resize(2);
	phi_struct[9].resize(1);

	int v1; int v2; int v3;
	for (int i = 0; i < 9; i++)
	{

		int index;
		index = (i > 2 && i < 9) ? 2 : 3;

		for (int j = 0; j < index; j++)
		{
			phi >> v1 >> v2 >> v3;
			phi_struct[i][j] = der_el(v1, v2, v3);
		}

	}
	for (int i = 2; i < 9; i++)
	{
		phi_struct[i][0].a = 13.5;
		phi_struct[i][1].a = -4.5;
	}
	for (int i = 0; i < 3; i++)
	{
		phi_struct[i][0].a = 4.5;
		phi_struct[i][1].a = -4.5;
		phi_struct[i][2].a = 1;
	}

	phi >> v1 >> v2 >> v3;
	phi_struct[9][0] = der_el(v1, v2, v3);
	phi_struct[9][0].a = 27;

	phi.close();
}

int fact(int N)
{
	int res = 1;
	for (int i = 2; i <= N; i++) res *= i;
	return res;
}

double Integral(vector<der_el> phii, vector<der_el> phij, double detD)
{
	int size_phii = phii.size;
	int size_phij = phij.size;
	double sum = 0;
	int v1, v2, v3;
	double coef = 0;

	for (int i = 0; i < size_phii; i++)
		for (int j = 0; j < size_phij; j++)
		{
			v1 = phii[i].v1 + phij[j].v1;
			v2 = phii[i].v2 + phij[j].v2;
			v3 = phii[i].v3 + phij[j].v3;
			coef = phii[i].a * phij[j].a;
			sum += coef * detD * (fact(v1) * fact(v2) * fact(v3)) / fact(v1 + v2 + v3 + 2);
		}
	return sum;
}

void G_Numbering()
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

void Portret()
{
	//связи связи( где для i базисной функции j не нулевая),

	//При прохождении по конечным элементам добавляем в каждый список
	//соответствующий ненулевой на этом элементе глобал баз ф-ции,
	//номера всех остальных ненулевых на этом элементе глобал базис ф-ций
	S.resize(global_number);// создаем множество где global_number - число глобал баз ф-ции 

	for (int t = 0; t < num_triangles; t++)////  проходим по конечным элементам (треугольники)
	{
		int index = 3;
		for (int i = 1; i < 10; i++) // проходим по каждой базисной функции кроме 1 
			for (int j = 0; j < i; j++, index += 2)// проходим по каждой базисной функции кроме 1 
			{

				auto a = triangle[t][i];// a - глобальный номер i , базисной функции на t-ом  конечном элементе 
				auto b = triangle[t][j];// b - вычисляем глобальный номер j , базисной функции на t-ом  конечном элементе 
				bool f = a < b; //проверяем какой глобальный номер больше 
				if (f)swap(a, b);//если f-true, меняем номера

				S[a].insert(b);//к a добавляем б
			}
	}
	//формируем портретматрицы (создаем массивы ig,jg)
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

void MxM(double detD, vector<vector<double>> local)
{

	for (int i = 0; i < 10; i++)
		for (int j = i; j < 10; j++)
		{
			local[i][j] = Integral(phi_struct[i], phi_struct[j], detD);
			//local[j][i] = local[i][j];
		}
}

void ReadDerivatives()
{
	ifstream in;
	in.open("derivatives.txt");
	for (int i = 0; i < 3; i++) derivatives[i].resize(3);
	for (int i = 3; i < 9; i++) derivatives[i].resize(4);
	derivatives[9].resize(3);

	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			in >> derivatives[i][j].gradNo
				>> derivatives[i][j].coeff
				>> derivatives[i][j].v1
				>> derivatives[i][j].v2
				>> derivatives[i][j].v3;
		}
	}

	for (int i = 3; i < 9; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			in >> derivatives[i][j].gradNo
				>> derivatives[i][j].coeff
				>> derivatives[i][j].v1
				>> derivatives[i][j].v2
				>> derivatives[i][j].v3;
		}
	}

	for (int j = 0; j < 3; j++)
	{
		in >> derivatives[9][j].gradNo
			>> derivatives[9][j].coeff
			>> derivatives[9][j].v1
			>> derivatives[9][j].v2
			>> derivatives[9][j].v3;
	}

	in.close();
}

void LocalPattern()
{
	localPattern.resize(10);
	for (int i = 0; i < 10; i++)
		localPattern[i].resize(10);

	for (int i = 0; i < 10; i++)
		for (int j = 0; j < 10; j++)
			for (auto x : derivatives[i])
				for (auto y : derivatives[j])
				{
					int v1 = x.v1 + y.v1;
					int v2 = x.v2 + y.v2;
					int v3 = x.v3 + y.v3;
					
					double coeff = x.coeff * y.coeff * fact(v1) * fact(v2) * fact(v3) / fact(v1 + v2 + v3 + 2);
					localPattern[i][j].push_back(local_el(x.gradNo, y.gradNo, coeff));
				}
}

void BuildGlobal()
{
	vector<vector<double>> global(global_number);
	for (int i = 0; i < global_number; i++)
		global[i].resize(global_number);

	for (int t = 0; t < num_triangles; t++)
	{
		auto tr = triangle[t];
		double D = detD(tr);
		alpha(tr);
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
				for (auto local : localPattern[i][j])
				{
					int i_global = tr[i];
					int j_global = tr[j];
					double scalGrad = grads[local.g1].first * grads[local.g2].first + grads[local.g1].second * grads[local.g2].second;
					global[i_global][j_global] += local.coeff * scalGrad * D;
				}
			}
		}
	}
}


int main()
{
	input();
	G_Numbering();
	Portret();




	return 0;
}