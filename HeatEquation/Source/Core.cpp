#include "Core.h"



//void FEM::InputBound(string bounds1File, string bounds2File)
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
//
//void FEM::Boundary2()
//{
//	vector<double> conditions;
//	conditions.resize(4);
//	edgeBasisValues.resize(4);
//
//	for (int edge = 0; edge < edgeCount2; edge++)
//	{
//		fill(conditions.begin(), conditions.end(), 0);
//		for (int i = 0; i < 4; i++)
//		{
//			double h = Distance(points[bound2[edge].v1], points[bound2[edge].v4]);
//			conditions[i] = h * edgeBasisValues[i] * thetaValue[bound2[edge].valueNo];
//		}
//
//		globalB[bound2[edge].v1] += conditions[0];
//		globalB[bound2[edge].v4] += conditions[1];
//		globalB[bound2[edge].v2] += conditions[2];
//		globalB[bound2[edge].v3] += conditions[3];
//	}
//}
//
//void FEM::Boundary1()
//{
//	for (auto edge : bound1)
//	{
//		globalMatrix.DI[edge.v1] = 1.0E+50;
//		globalB[edge.v1] = 1.0E+50 * uValue[edge.valueNo];
//		globalMatrix.DI[edge.v2] = 1.0E+50;
//		globalB[edge.v2] = 1.0E+50 * uValue[edge.valueNo];
//		globalMatrix.DI[edge.v3] = 1.0E+50;
//		globalB[edge.v3] = 1.0E+50 * uValue[edge.valueNo];
//		globalMatrix.DI[edge.v4] = 1.0E+50;
//		globalB[edge.v4] = 1.0E+50 * uValue[edge.valueNo];
//	}
//}
