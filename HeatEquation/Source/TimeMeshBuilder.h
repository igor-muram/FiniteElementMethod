#pragma once

#include <iostream>
#include <vector>
#include "Interval.h"

using TimeIterator = std::vector<double>::iterator;

class TimeMeshBuilder
{
public:
	TimeMeshBuilder(std::vector<Interval>& intervals)
	{
		for (auto interval : intervals)
		{
			double begin = interval.begin;
			double end = interval.end;
			int n = interval.n;
			double q = interval.q;

			double h = (q == 1.0) ? ((end - begin) / n) : ((end - begin) * (1.0 - q) / (1.0 - pow(q, n)));

			t.push_back(begin);
			for (int i = 1; i < n; i++)
			{
				t.push_back(t.back() + h);
				h *= q;
			}
		}

		t.push_back(intervals.back().end);
	}

	const std::vector<double>& getMesh() { return t; }

	TimeIterator Begin() { return t.begin(); }
	TimeIterator End() { return t.end(); }

	int Size() { return t.size(); }

private:
	std::vector<double> t;
};