// ChiaSig: detection of significant interactions in ChIA-PET data using Non-central hypergeometric distribution (NCHG)
// This file is subject to the terms and conditions defined in file 'LICENSE', which is part of this source code package.
// Copyright (2014) Jonas Paulsen


#ifndef GUARD_CCCMatrix
#define GUARD_CCCMatrix

#include "Segment.h"
#include <vector>
#include <set>
#include <map>
#include <utility>
#include <limits>
#include <iostream>

#include "../spline/spline.h" // Spline-smoothing

const int MAXDELTA = std::numeric_limits<int>::max();

typedef int32_t SegmentKey;
enum IndexType {NONINDEX, INTRAINDEX, ROWINDEX, COLINDEX};

#define CACHE_FLAG_N       0x00000001
#define CACHE_FLAG_EXPECTN 0x00000002



class CCCMatrixInteraction {
public:
	CCCMatrixInteraction();
	CCCMatrixInteraction(std::set<Segment>&, std::set<Segment>&, Interaction defaultvalue, bool, std::string& chr);
	~CCCMatrixInteraction();
	void InteractionLoaded(int minDelta/*=8000*/, int maxDelta/*=MAXDELTA*/);

	uint32_t getCountElement(SegmentKey, SegmentKey);
	void setElement(Segment&, Segment&, Interaction&);
	void getRowSums(std::vector<uint32_t > &);
	
	// TODO: decide whether these are external aux function or members
	void getStatisticsPerDelta(std::map<int, DeltaStatistics >&, QuantileMapper*, int minDelta=8000, int maxDelta=MAXDELTA, bool masking=false);

	int getNumberOfPositives(double alpha);
	int getNumberOfPositives(double alpha, int cutoff);
	unsigned long maskByAlpha(double alpha, int cutoff, bool dump=false);
	void calculatePvalues(int minDelta=8000, int maxDelta=MAXDELTA);
	void printPositives(std::string&, double alpha, int cutoff, bool printNij, bool printAll=false);
	// TODO: END - decide whether these are external aux function or members
	
	unsigned long getDeltasCount(int minDelta /* = 0 */, int maxDelta /*=MAXDELTA*/, bool masking /* =false */);
	unsigned long getDeltasCount(DeltaIndicatorCounter*, int minDelta /* = 0 */, int maxDelta /*=MAXDELTA*/, bool masking /* =false */);
	uint32_t getN(int minDelta=0, int maxDelta=MAXDELTA) {
		return cachedN;
	}
	bool isIntra;
	
	// TODO: we are keeping the expectation values
	void setDeltaToExpectationMapper(spline& smoothDelta2ExpMapper, int minDelta=8000, int maxDelta=MAXDELTA);
	float getExpectN(int minDelta=0, int maxDelta=MAXDELTA) {
		return cachedExpectN;
	}
	float getExpectElement(SegmentKey rowKey, SegmentKey colKey) {
		// TODO:
		int delta = abs(orderedSegmentMidPoints[colKey]-orderedSegmentMidPoints[rowKey]);
		if(delta>=minExpDelta and delta<=maxExpDelta) {
			return pDelta2ExpSpline->cachedAt(delta);
		} else {
			return 0.0f;
		}
	}
	void printExpectationMatrix(std::ostream&, bool hideMask=false);
	void getExpectRowSums(std::vector<float > &);
	// TODO: END - we are keeping the expectation values
	
	// TODO: for FDR
	double estimateFalsePositives(double alpha, int cutoff, int minDelta=8000, int maxDelta=MAXDELTA);
	uint32_t getSegmentCount();
	// TODO: END - for FDR
		
	// TODO: assess caching effect
	void cacheRowSums();
	// TODO: END - assess caching effect
	
	unsigned long getInteractionCount();
	void setThreads(int n);
	
	// for deltas managed memory
	int getWidestInteractionSpan();
	// END - for deltas managed memory
	
private:
	Interaction defaultVal;
	std::map<Segment, SegmentKey > segmentToKey;
	std::vector<SegmentMin> orderedSegments;
	int* orderedSegmentMidPoints;
	std::vector<std::map<SegmentKey, Interaction > > matrix;
	
	// TODO: for reporting
	unsigned long m_TotalInteractions;
	
	// TODO: assess caching effect
	std::vector<float > m_getExp;
	std::vector<uint32_t > m_getMarginal;
	// TODO: END - assess caching effect

	
	// TODO: we are keeping the expectation values
	spline* pDelta2ExpSpline;
	int minExpDelta;
	int maxExpDelta;
	// TODO: END - we are keeping the expectation values
	
	int nThreads;
	
	// TODO: more caching flags
	uint32_t cacheFlag;
	uint32_t cachedN;
	void calculate_cache_N(int minDelta=0, int maxDelta=MAXDELTA);
	float cachedExpectN;
	void calculate_cache_ExpectN(int minDelta=0, int maxDelta=MAXDELTA);
	
	// TODO: for easier debugging
	std::string m_chr;
};

#endif
