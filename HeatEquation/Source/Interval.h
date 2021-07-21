#pragma once

#include <iostream>

struct Interval
{
	double begin = 0.0, end = 0.0;
	double q = 0.0;
	int n = 0;
};

inline std::istream& operator>>(std::istream& stream, Interval& interval)
{
	stream >> interval.begin >> interval.end >> interval.n >> interval.q;
	return stream;
}