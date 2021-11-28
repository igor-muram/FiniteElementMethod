#pragma once

#include <vector>
#include <set>

#include "FEMInfo.h"
#include "Matrix.h"

using namespace std;

class PortraitBuilder
{
public:
	PortraitBuilder(int nodeCount, vector<FiniteElement>& elements) : nodeCount(nodeCount)
	{
		connections.resize(nodeCount);
		BuildConnections(elements);
	}

	void Build(Matrix& A)
	{
		A.N = nodeCount;
		A.IA.resize(nodeCount + 1);
		A.IA[0] = A.IA[1] = 0;

		for (int i = 2; i <= nodeCount; i++)
		{
			int col = A.IA[i - 1];
			A.IA[i] = col + connections[i - 1].size();
		}

		A.JA.resize(A.IA[nodeCount]);

		for (int i = 1, k = 0; i < nodeCount; i++)
		{
			for (int j : connections[i])
			{
				A.JA[k] = j;
				k++;
			}
		}

		A.DI.resize(nodeCount);
		A.AL.resize(A.IA[nodeCount]);
	}

private:
	int nodeCount, JASize;
	vector<set<int>> connections;

	void BuildConnections(vector<FiniteElement>& elements)
	{
		for (auto& t : elements)
			for (int i = 1; i < basisSize; i++)
				for (int j = 0; j < i; j++)
				{
					auto a = t.verts[i];
					auto b = t.verts[j];
					if (a < b) swap(a, b);

					connections[a].insert(b);
				}
	}
};