#pragma once

#include <vector>
#include <map>
#include "FEMInfo.h"

using namespace std;

class MeshBuilder
{
public:
	MeshBuilder(int nodeCount): nodeCount(nodeCount)
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

public:
	int nodeCount = 0;
	vector<map<int, int>> edgeMatrix;
};