// ChiaSig: detection of significant interactions in ChIA-PET data using Non-central hypergeometric distribution (NCHG)
// This file is subject to the terms and conditions defined in file 'LICENSE', which is part of this source code package.
// Copyright (2014) Jonas Paulsen

#ifndef GUARD_Segment
#define GUARD_Segment

#include <string>
#include <map> // QuantileMapper
#include <vector> // QuantileMapper
#include <stdlib.h> //abs
#include <stdint.h> //uint16_t

typedef uint16_t ChromosomeIndexType;
#include "../spline/spline.h"

struct Segment {
	Segment();
	Segment(ChromosomeIndexType c, int s, int e);
	
	ChromosomeIndexType chr;
	int start;
	int end;
	int pos;
	
	bool operator < (const Segment& str) const {
		return (chr == str.chr ? pos < str.pos : chr < str.chr);
	}
	bool operator == (const Segment& str) const {
		return (chr == str.chr and start == str.start and end == str.end);
	}
	
	int operator - (const Segment& str) const {
		return (chr == str.chr ? abs(pos - str.pos) : -1);
	}
};

struct SegmentMin {
	SegmentMin();
	SegmentMin(ChromosomeIndexType c, int s, int e);
	SegmentMin(const Segment &);
	
	ChromosomeIndexType chr;
	int start;
	int end;
};

struct Interaction {
	Interaction();
	Interaction(uint32_t c, bool m, double p);
	uint32_t mask:1, count:31;
	double pvalue;
};

struct DeltaStatistics {
	DeltaStatistics();
	DeltaStatistics(uint32_t s, uint32_t c);
	float getMean() { return (sum*1.0f/count); }
	uint32_t sum;
	uint32_t count;
};


struct QuantileMapper {
	QuantileMapper() { }
	void computeDeltaCountQuantile(unsigned long, DeltaIndicatorCounter*, unsigned int, bool masking=false);
	std::vector<int>& getQuantileDeltas();
	int getDeltaQuantile(int);
	std::vector<int> m_deltas;
	std::vector<int> m_quantiles;
};

#endif
