// ChiaSig: detection of significant interactions in ChIA-PET data using Non-central hypergeometric distribution (NCHG)
// This file is subject to the terms and conditions defined in file 'LICENSE', which is part of this source code package.
// Copyright (2014) Jonas Paulsen


#include "Segment.h"

#include <iostream>
#include <fstream>
#include <algorithm>
#include <numeric>

#include <assert.h>

using namespace std;

extern "C" {
	double cputime();
	double realtime();
}

Segment::Segment()
:chr(0),start(0),end(0),pos(0)
{
}

Segment::Segment(ChromosomeIndexType c, int s, int e)
:chr(c), start(s), end(e), pos((s+e)/2)
{
}

SegmentMin::SegmentMin()
:chr(0),start(0),end(0)
{
}

SegmentMin::SegmentMin(ChromosomeIndexType c, int s, int e)
:chr(c), start(s), end(e)
{
}

SegmentMin::SegmentMin(const Segment &a)
:chr(a.chr), start(a.start), end(a.end)
{
}


Interaction::Interaction()
:count(0), mask(0), pvalue(1)
{
}

Interaction::Interaction(uint32_t c, bool m, double p)
:count(c), mask(m ? 1 : 0), pvalue(p)
{
}

DeltaStatistics::DeltaStatistics()
: sum(0), count(0)
{
}

DeltaStatistics::DeltaStatistics(uint32_t s, uint32_t c)
: sum(s), count(c)
{
}

// WCH: this version keeps the counting of the deltas instead of vector of deltas
void QuantileMapper::computeDeltaCountQuantile(unsigned long vecSize, DeltaIndicatorCounter* vec, unsigned int nQuant, bool masking)
{
	m_deltas.clear(); m_deltas.reserve(nQuant+1);
	m_quantiles.clear(); m_quantiles.reserve(nQuant+1);
	
	// Function that computes quantiles of size quantSize and returns
	// a map from each integer element to its quantile.
	// Obs: If a given for delta is found in two adjacent divisions,
	// then the smallest delta is used.
	
	double t_diff;
	double t_start;
	
	vector<int>::size_type quantSize = vecSize/nQuant;
	fprintf(stderr, "[M::%s:%d] vecSize=%lu, nQuant=%d, quantSize=%lu\n", __func__, __LINE__, vecSize, nQuant, quantSize);
	assert(vecSize > quantSize);
	
	t_start = realtime();
	vector<int>::size_type nTotalQuantSize = quantSize*nQuant;
	int nProcessedWRTdata = 0;
	if (nTotalQuantSize>vecSize) nProcessedWRTdata = 1;
	else if (nTotalQuantSize<vecSize) nProcessedWRTdata = -1;
	// 0==nProcessedWRTdata : exact bin size matched
	// nProcessedWRTdata < 0 : i.e. quantSize*nQuant < vecSize
	//                         we have some spill over, which we will just used the last chunk media
	// nProcessedWRTdata > 0 : i.e. quantSize*nQuant > vecSize
	//                         last chunk will have insufficient element in bin
	int nEffQuant = (nProcessedWRTdata > 0) ? nQuant - 1 : nQuant;
	int median = 0;
	vector<int>::size_type i = 0, j=0;

	// loop thru' once just for the indexes so that we can remember the values that we need
	// then, we look up values from these needed position for the actual computation
	// we should also keep the 1st and last deltas
	vector<unsigned long> indices; //indices.push_back(0);
	for (int nQuantIndex=0; nQuantIndex < nEffQuant; ++nQuantIndex) {
		indices.push_back(i); // before debugging dump
		vector<int>::size_type midpos=i+(quantSize)/2;
		if (quantSize % 2) {
			indices.push_back(midpos);
		} else {
			indices.push_back(midpos-1);
			indices.push_back(midpos);
		}
		
		j = i + (quantSize-1);
		indices.push_back(j);
		
		i+=quantSize;
	}
	// Last remaindig chunk of values:
	if (0!=nProcessedWRTdata) {
		indices.push_back(i); // before debugging dump
		vector<int>::size_type tmp = quantSize;
		quantSize=(vecSize-i); // quantSize of final chunk
		j = i + (quantSize-1);
		
		if (0>nProcessedWRTdata) {
			// handle the spill over
			// WCH: IMPORTANT: we are mimicking the original implementation!!!
			//      Since the map is "full", set to the last value
			
			// optimized implementation
			indices.push_back(j);
		} else {
			// handle insufficient data in last bin
			vector<int>::size_type midpos=i+(quantSize)/2;
			if (quantSize % 2) {
				indices.push_back(midpos);
			} else {
				indices.push_back(midpos-1);
				indices.push_back(midpos);
			}
			
			indices.push_back(j);
		}
		quantSize = tmp;
	}
	indices.push_back(vecSize-1);
	
	// we now go thru' the deltasCount vector to cache the values at the required indices
	sort(indices.begin(), indices.end());
	map<unsigned long, int > vecMap;
	unsigned long nCurrCount = 0;
	unsigned long indicesIdx = 0;
	for(int i=0; nCurrCount<vecSize && indicesIdx<indices.size(); ++i) {
		if (DELTA_COUNTER(vec[i])>0) {
			nCurrCount += DELTA_COUNTER(vec[i]);
			// (nCurrCount-1) for getting the last index, as index is 0-based
			while ((nCurrCount-1)>=indices[indicesIdx]) {
				// passed a band
				vecMap[indices[indicesIdx]] = i;
				indicesIdx++;
				if (indicesIdx>=indices.size()) {
					break;
				}
			}
		}
	}
	// END - hacked
	
	// re-initialized as the previous section has tamed the variable
	i = 0;
	// let's build the look up table
	for (int nQuantIndex=0; nQuantIndex < nEffQuant; ++nQuantIndex) {
		vector<int>::size_type midpos=i+(quantSize)/2;
		median = (quantSize % 2) ? vecMap[midpos] : (vecMap[midpos-1]+vecMap[midpos])/2;
		
		j = i + (quantSize-1);
		if (0==m_deltas.size()) {
			m_deltas.push_back(vecMap[j]);
			m_quantiles.push_back(median);
		} else if (m_deltas.back()!=vecMap[j]) {
			m_deltas.push_back(vecMap[j]);
			m_quantiles.push_back(median);
		}
		
		i+=quantSize;
	}
	
	int lastval=median;
	// Last remaindig chunk of values:
	if (0!=nProcessedWRTdata) {
		quantSize=(vecSize-i); // quantSize of final chunk
		j = i + (quantSize-1);
		
		if (0>nProcessedWRTdata) {
			// handle the spill over
			// WCH: IMPORTANT: we are mimicking the original implementation!!!
			//      Since the map is "full", set to the last value
			
			// optimized implementation
			if (0==m_deltas.size()) {
				m_deltas.push_back(vecMap[j]);
				m_quantiles.push_back(lastval);
			} else {
				m_deltas.back() = vecMap[j];
			}
			
		} else {
			// handle insufficient data in last bin
			vector<int>::size_type midpos=i+(quantSize)/2;
			median = (quantSize % 2) ? vecMap[midpos] : (vecMap[midpos-1]+vecMap[midpos])/2;
			
			if (0==m_deltas.size()) {
				m_deltas.push_back(vecMap[j]);
				m_quantiles.push_back(median);
			} else if (m_deltas.back()!=vecMap[j]) {
				m_deltas.push_back(vecMap[j]);
				m_quantiles.push_back(median);
			}
			
		}
	}
	
	t_diff = realtime() - t_start;
	fprintf(stderr, "[M::%s:%d] setup took %.2f min..\n", __func__, __LINE__, t_diff/60.0);
}

vector<int>& QuantileMapper::getQuantileDeltas()
{
	return m_deltas;
}

int QuantileMapper::getDeltaQuantile(int delta)
{
	// check look up table
	vector<int>::iterator it = lower_bound(m_deltas.begin(),m_deltas.end(),delta);
	int idx = std::max( int(it-m_deltas.begin()), 0);
	return m_quantiles[idx];
}
