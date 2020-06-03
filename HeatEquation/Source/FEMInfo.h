#pragma once

#include <vector>
#include <functional>

using namespace std;

struct FiniteElement
{
	std::vector<int> verts;
	int materialNo;
};

const int basisSize = 10;

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

vector<vector<DerComp>> ders = {
	{ { 0, 13.5, 2, 0, 0 },			{ 0, -9, 1, 0, 0 },			{ 0, 1, 0, 0, 0 } },
	{ { 1, 13.5, 0, 2, 0 },			{ 1, -9, 0, 1, 0 },			{ 1, 1, 0, 0, 0 } },
	{ { 2, 13.5, 0, 0, 2 },			{ 2, -9, 0, 0, 1 },			{ 2, 1, 0, 0, 0 } },
	{ { 0, 27, 1, 1, 0 },			{ 0, -4.5, 0, 1, 0 },		{ 1, 13.5, 2, 0, 0 },		{ 1, -4.5, 1, 0, 0 } }, // 4
	{ { 0, 13.5, 0, 2, 0 },			{ 0, -4.5, 0, 1, 0 },		{ 1, 27, 1, 1, 0 },			{ 1, -4.5, 1, 0, 0 } }, // 5

	{ { 0, 27, 1, 0, 1 },			{ 0, -4.5, 0, 0, 1 },		{ 2, 13.5, 2, 0, 0 },		{ 2, -4.5, 1, 0, 0 } }, // 9
	{ { 0, 13.5, 0, 0, 2 },			{ 0, -4.5, 0, 0, 1 },		{ 2, 27, 1, 0, 1 },			{ 2, -4.5, 1, 0, 0 } }, // 8

	{ { 1, 27, 0, 1, 1 },			{ 1, -4.5, 0, 0, 1 },		{ 2, 13.5, 0, 2, 0 },		{ 2, -4.5, 0, 1, 0 } }, // 6
	{ { 1, 13.5, 0, 0, 2 },			{ 1, -4.5, 0, 0, 1 },		{ 2, 27, 0, 1, 1 },			{ 2, -4.5, 0, 1, 0 } }, // 7
	{ { 0, 27, 0, 1, 1 },			{ 1, 27, 1, 0, 1 },			{ 2, 27, 1, 1, 0 } }
};

vector<vector<PsiComp>> basis = {
	{ { 4.5, 3, 0, 0 },         { -4.5, 2, 0, 0 },       { 1, 1, 0, 0 } },
	{ { 4.5, 0, 3, 0 },         { -4.5, 0, 2, 0 },       { 1, 0, 1, 0 } },
	{ { 4.5, 0, 0, 3 },         { -4.5, 0, 0, 2 },       { 1, 0, 0, 1 } },
	{ { 13.5, 2, 1, 0 },        { -4.5, 1, 1, 0 } },
	{ { 13.5, 1, 2, 0 },        { -4.5, 1, 1, 0 } },

	{ { 13.5, 2, 0, 1 },        { -4.5, 1, 0, 1 } },
	{ { 13.5, 1, 0, 2 },        { -4.5, 1, 0, 1 } },

	{ { 13.5, 0, 2, 1 },        { -4.5, 0, 1, 1 } },
	{ { 13.5, 0, 1, 2 },	    { -4.5, 0, 1, 1 } },
	{ { 27, 1, 1, 1 } }
};

vector<function<double(double, double, double)>> f = {
	[](double x, double y, double t) { return 2 * t; }
};

vector<function<double(double, double, double)>> uValue = {
	[](double x, double y, double t) { return t * t; }
};
vector<double> thetaValue = { 0, 0, 0 };
vector<double> edgeBasisValues = { 0.125, 0.125, 0.375, 0.375 };

double u(double x, double y, double t)
{
	return t * t;
}