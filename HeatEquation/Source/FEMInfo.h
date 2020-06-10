#pragma once

#include <vector>
#include <functional>

using namespace std;

struct FiniteElement
{
	vector<int> verts;
	int materialNo = 0;
};

const int basisSize = 10;

struct DerComp
{
	int gradNo = 0;
	double coeff = 0.0;
	int v1 = 0, v2 = 0, v3 = 0;
};

struct PsiComp
{
	double coeff = 0.0;
	int v1 = 0, v2 = 0, v3 = 0;
};

struct LocalComp
{
	LocalComp(int grad1, int grad2, double coeff) : grad1(grad1), grad2(grad2), coeff(coeff) {}
	int grad1 = 0, grad2 = 0;
	double coeff = 0.0;
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
	{ { 13.5, 0, 1, 2 },		{ -4.5, 0, 1, 1 } },
	{ { 27, 1, 1, 1 } }
};

vector<function<double(double, double, double)>> uValue = {
	[](double x, double y, double t) { return t * t; }
};

vector<function<double(double, double, double)>> thetaValue = {
	[](double x, double y, double t) { return 0; },
	[](double x, double y, double t) { return 0; }
};

double u(double x, double y, double t)
{
	return t * t;
}