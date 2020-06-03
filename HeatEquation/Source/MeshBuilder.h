#pragma once

#include <vector>
#include <map>
#include "FEMInfo.h"
#include "Boundary.h"

using namespace std;

class MeshBuilder
{
public:
	MeshBuilder(int nodeCount) : nodeCount(nodeCount)
	{
		edgeMatrix.resize(nodeCount);
	}

	void Build(std::vector<FiniteElement>& elements)
	{
		for (auto& e : elements)
		{
			int index = 3;
			for (int i = 0; i < 3; i++)
				for (int j = i + 1; j < 3; j++, index += 2)
				{
					int a = e.verts[i];
					int b = e.verts[j];
					bool f = a > b;
					if (f) std::swap(a, b);

					if (edgeMatrix[a][b] == 0)
					{
						e.verts[index] = nodeCount + (f ? 1 : 0);
						e.verts[index + 1] = nodeCount + (f ? 0 : 1);
						edgeMatrix[a][b] = nodeCount;
						nodeCount += 2;
					}
					else
					{
						int a1 = edgeMatrix[a][b];
						e.verts[index] = a1 + (f ? 1 : 0);
						e.verts[index + 1] = a1 + (f ? 0 : 1);
					}
				}
			e.verts[index] = nodeCount;
			nodeCount++;
		}
	}

	void BuildBoundary(vector<Edge>& edges)
	{
		for (auto& edge : edges)
		{
			auto a = edge.v1;
			auto b = edge.v4;
			bool f = a > b;
			if (f) swap(a, b);

			edge.v2 = edgeMatrix[a][b] + (f ? 0 : 1);
			edge.v3 = edgeMatrix[a][b] + (f ? 1 : 0);
		}
	}

public:
	int nodeCount = 0;
	vector<map<int, int>> edgeMatrix;
};