// ChiaSig: detection of significant interactions in ChIA-PET data using Non-central hypergeometric distribution (NCHG)
// This file is subject to the terms and conditions defined in file 'LICENSE', which is part of this source code package.
// Copyright (2014) Jonas Paulsen


#include "CCCMatrix.h"
#include <assert.h>
#include <iostream>
#include <fstream>  // for logging
#include <numeric>  // for accumulate
#include "../stocc/stocc.h" // Non-central hypergeometric and Poisson
#include <string.h> // for memcpy

void kt_for(int n_threads, void (*func)(void*,int,int), void *data, int n);

using namespace std;

const float fSmallValue = 1E-20;
const double dOmegaPrecision = 1E-50;
const double dOmegaTablePrecision = 1E-53; //dOmegaPrecision * 0.001

const double dMinOmega=1E-01;

#define ORIGINAL_CODE
#ifdef ORIGINAL_CODE
// original

//double
//   getQuantilePvalNCHG(int nij, double p, int ni, int nj, int n, double odds)
double getCumulativeNCHG(int nij, int ni, int nj, int N, double odds, double precision /*=1E-20*/, bool lower_tail /*=true*/) {
	// Some assertions here
	int n=nj;
	int m1=ni;
	//int m2=N-ni;
	
	int xmin = m1 + n - N;  if (xmin < 0) xmin = 0;  // Minimum x
	int xmax = n;  if (xmax > m1) xmax = m1;         // Maximum x
	
	//int BufferLength;
	int x1, x2, x;
	double sum, p;
	double* buff = 0;
	
	CFishersNCHypergeometric xhyp(nj, ni, N, odds, precision);
	int BufferLength = (int)xhyp.MakeTable(buff, 0, &x1, &x2, precision * 0.001);
	double buffer[BufferLength];
	double factor = 1. / xhyp.MakeTable(buffer, BufferLength, &x1, &x2, precision * 0.001); // make table of probability values
	int xmean = (int)(xhyp.mean() + 0.5);           // Round mean
	
	assert(xmean >= x1);
	assert(xmean <= x2);
	
	for (x = x1, sum = 0; x <= xmean; x++) sum = buffer[x-x1] += sum; // Make left tail of table cumulative
	for (x = x2, sum = 0; x > xmean; x--) sum = buffer[x-x1] += sum;  // Make right tail of table cumulative from the right
	
	
	x = nij;                       // Input x value
	if (x <= xmean) { // Left tail
		if (x < x1) {
			p = 0.;                    // Outside table
		}
		else {
			p = buffer[x-x1] * factor; // Probability from table
		}
		if (!lower_tail) p = 1. - p;  // Invert if right tail
	}
	else {      // Right tail
		if (x >= x2) {
			p = 0.;                    // Outside table
		}
		else {
			p = buffer[x-x1+1] * factor; // Probability from table
		}
		if (lower_tail) p = 1. - p;   // Invert if left tail
	}
	return p;
}

//double
//  getPvalNCHGNew(int nij, int ni, int nj, int n, double odds)
double getPvalNCHG(int nij, int ni, int nj, int n, double odds, double precision /*=1E-20*/, double minOmega /*=0.1*/) {
	if(nij>0) {
		return getCumulativeNCHG(nij-1, ni, nj, n, max(odds, minOmega), precision, false);
	}
	else {
		return 1;
	}
}

//int
//  getQuantileNCHG(double p, int ni, int nj, int n, double odds)
int getQuantileNCHG(double p, int ni, int nj, int n, double odds, double precision, bool lower_tail /*=true*/) {
	assert(p >= 0 and p <=1);
	//int BufferLength;
	int x=0;
	int x1=0;
	int x2=0;
	double* buff = 0;
	double sum;
	unsigned int a, b, c; // Used in binary search
	
	CFishersNCHypergeometric xhyp(nj, ni, n, odds, precision);
	int BufferLength = (int)xhyp.MakeTable(buff, 0, &x1, &x2, precision * 0.001);
	if(BufferLength == 0) { return -1; } // Buffer cannot be allocated
	double buffer[BufferLength];
	double factor = xhyp.MakeTable(buffer, BufferLength, &x1, &x2, precision * 0.001); // make table of probability values
	
	for (x = x1, sum = 0; x <= x2; x++) sum = buffer[x-x1] += sum; // Make table cumulative
	
	
	if (!lower_tail) p = 1. - p; // Invert if right tail
	p *= factor; // Table is scaled by factor
	
	// Binary search in table:
	a = 0; b = x2 - x1 + 1;
	while (a < b) {
		c = (a + b) / 2;
		if (p <= buffer[c]) {
			b = c;
		}
		else {
			a = c + 1;
		}
	}
	x = x1 + a;
	if (x > x2) x = x2;           // Prevent values > xmax that occur because of small imprecisions
	return x;
}

#else
// WCH : OPTIMZED
double getPvalNCHGNew(int nij, int ni, int nj, int n, double odds)
{
	nij--;
	odds = max(odds, dMinOmega);
	
	CFishersNCHypergeometric xhyp(nj, ni, n, odds, dOmegaPrecision);
	int BufferLength = xhyp.getTableLength();
	double buffer[BufferLength];
	int x1, x2;
	double factor = 1. / xhyp.MakeTable(buffer, BufferLength, &x1, &x2, dOmegaTablePrecision); // make table of probability values
	int xmean = (int)(xhyp.mean() + 0.5);           // Round mean
	
	assert(xmean >= x1);
	assert(xmean <= x2);
	
	double dval;
	if (nij <= xmean) { // Left tail
		if (nij < x1) {
			dval = 1.;
		} else {
			double sum; int x;
			for (x = x1, sum = 0; x <= nij; x++) sum = buffer[x-x1] += sum; // Make left tail of table cumulative
			dval = 1. - buffer[nij-x1] * factor;
		}
	} else {      // Right tail
		if (nij >= x2) {
			dval = 0.;
		} else {
			double sum; int x;
			for (x = x2, sum = 0; x > nij; x--) sum = buffer[x-x1] += sum;  // Make right tail of table cumulative from the right
			dval = buffer[nij-x1+1] * factor;
		}
	}
	return dval;
}
// END - WCH : OPTIMZED

// collapse of getQuantileNCHG()+getPvalNCHGNew()
double getQuantilePvalNCHG(int nij, double p, int ni, int nj, int n, double odds)
{
	
	assert(p >= 0 and p <=1);
	assert(odds>=dMinOmega);
	
	CFishersNCHypergeometric xhyp(nj, ni, n, odds, dOmegaPrecision);
	int BufferLength = xhyp.getTableLength();
	if(BufferLength == 0) { return -1; } // Buffer cannot be allocated
	double buffer[BufferLength];
	int x1,x2;
	double factor = xhyp.MakeTable(buffer, BufferLength, &x1, &x2, dOmegaTablePrecision); // make table of probability values
	double secbuffer[BufferLength];
	memcpy(secbuffer, buffer, BufferLength*sizeof(double));
	
	int x;
	double sum;
	for (x = x1, sum = 0; x <= x2; x++) sum = buffer[x-x1] += sum; // Make table cumulative
	
	p *= factor; // Table is scaled by factor
	
	// Binary search in table:
	unsigned int a, b, c; // Used in binary search
	a = 0; b = x2 - x1 + 1;
	while (a < b) {
		c = (a + b) / 2;
		if (p <= buffer[c]) {
			b = c;
		}
		else {
			a = c + 1;
		}
	}
	x = x1 + a;
	if (x > x2) x = x2;           // Prevent values > xmax that occur because of small imprecisions
	
	nij = max(nij, 1+x);
	nij--; //nij-1
	
	factor = 1. / factor;
	int xmean = (int)(xhyp.mean() + 0.5);           // Round mean
	
	assert(xmean >= x1);
	assert(xmean <= x2);
	
	// NOTE: no need to re-assign, just use nij directly
	double dval;
	if (nij <= xmean) { // Left tail
		if (nij < x1) {
			dval =  1.;
		} else {
			for (x = x1, sum = 0; x <= nij; x++) sum = secbuffer[x-x1] += sum; // Make left tail of table cumulative
			dval = 1. - secbuffer[nij-x1] * factor;
		}
	} else {      // Right tail
		if (nij >= x2) {
			dval = 0.;
		} else {
			for (x = x2, sum = 0; x > nij; x--) sum = secbuffer[x-x1] += sum;  // Make right tail of table cumulative from the right
			dval = secbuffer[nij-x1+1] * factor;
		}
	}
	return dval;
}

// original version
int getQuantileNCHG(double p, int ni, int nj, int n, double odds)
{
	assert(p >= 0 and p <=1);
	
	CFishersNCHypergeometric xhyp(nj, ni, n, odds, dOmegaPrecision);
	int BufferLength = xhyp.getTableLength();
	if(BufferLength == 0) { return -1; } // Buffer cannot be allocated
	double buffer[BufferLength];
	int x1,x2;
	double factor = xhyp.MakeTable(buffer, BufferLength, &x1, &x2, dOmegaTablePrecision); // make table of probability values
	
	int x;
	double sum;
	for (x = x1, sum = 0; x <= x2; x++) sum = buffer[x-x1] += sum; // Make table cumulative
	
	p *= factor; // Table is scaled by factor
	
	// Binary search in table:
	unsigned int a, b, c; // Used in binary search
	a = 0; b = x2 - x1 + 1;
	while (a < b) {
		c = (a + b) / 2;
		if (p <= buffer[c]) {
			b = c;
		}
		else {
			a = c + 1;
		}
	}
	x = x1 + a;
	if (x > x2) x = x2;           // Prevent values > xmax that occur because of small imprecisions
	
	return x;
}
#endif

//------------------------------------------------------------------------------


CCCMatrixInteraction::CCCMatrixInteraction()
:isIntra(true),defaultVal(Interaction()),m_TotalInteractions(0),nThreads(1)
,orderedSegmentMidPoints(NULL)
,cacheFlag(0),cachedN(0),cachedExpectN(0.0)
{
}

CCCMatrixInteraction::CCCMatrixInteraction(set<Segment> &rows, set<Segment> &cols, Interaction defaultvalue, bool isintra, string& chr)
:isIntra(isintra),defaultVal(defaultvalue),m_TotalInteractions(0),nThreads(1)
,orderedSegmentMidPoints(NULL)
,cacheFlag(0),cachedN(0),cachedExpectN(0.0)
,m_chr(chr)
{
	// Sanity checks:
	assert(rows.size()>0);
	assert(cols.size()>0);
	if (isintra) {
		assert(rows.size()==cols.size());
	}
	else {
		// Obs: This could be made more sophisticated.
		// None of the elements in rows and cols should be allowed to be the same!
		assert(rows.size()!=cols.size());
	}
	
	assert(isintra); // WCH: we will have to handle trans-interaction differently (future work)
	// create the look up
	SegmentKey segmentCounter=0;
	for (set<Segment>::iterator it = rows.begin(); it != rows.end(); ++it) {
		segmentToKey[*it] = segmentCounter;
		segmentCounter++;
	}
	// store as row-major
	matrix.resize(segmentCounter);
}

CCCMatrixInteraction::~CCCMatrixInteraction()
{
	if (orderedSegmentMidPoints) {
		free(orderedSegmentMidPoints);
		orderedSegmentMidPoints = NULL;
	}
}

void CCCMatrixInteraction::InteractionLoaded(int minDelta/*=8000*/, int maxDelta/*=MAXDELTA*/)
{
	// all the contacts have been loaded
	// we will switch over to crunching of data only
	if (0==orderedSegments.size()) {
		orderedSegments.resize(segmentToKey.size());
		orderedSegmentMidPoints = (int *) calloc(segmentToKey.size(), sizeof(int));
		if (!orderedSegmentMidPoints) {
			fprintf(stderr, "[E::%s:%d] Fail to allocate %lu unit of midpoint.\n", __func__, __LINE__, segmentToKey.size());
			exit(-1);
		}
		for (map<Segment, SegmentKey >::iterator itRow = segmentToKey.begin(); itRow != segmentToKey.end(); ++itRow) {
			orderedSegments[itRow->second] = SegmentMin(itRow->first);
			orderedSegmentMidPoints[itRow->second] = itRow->first.pos;
		}
		// we no longer needs the keys
		segmentToKey.clear();
	}
}

// NO optimization needed
uint32_t CCCMatrixInteraction::getCountElement(SegmentKey rowKey, SegmentKey colKey) {
	map<SegmentKey, Interaction >::iterator it = matrix[rowKey].find(colKey);
	if (it==matrix[rowKey].end()) {
		return defaultVal.count;
	} else {
		return it->second.count;
	}
}

// NO optimization needed
// WCH: this is the only version to set the element in a matrix by "Segment" rather than "SegmentKey"
void CCCMatrixInteraction::setElement(Segment& row, Segment& col, Interaction& val) {
	SegmentKey rowKey = segmentToKey.at(row);
	SegmentKey colKey = segmentToKey.at(col);
	
	assert(rowKey!=colKey);
	if (rowKey > colKey) {
		SegmentKey t = colKey; colKey = rowKey; rowKey = t;
	}
	
	m_TotalInteractions++;
	map<SegmentKey, Interaction >::iterator it = matrix[rowKey].find(colKey);
	if (it==matrix[rowKey].end()) {
		matrix[rowKey][colKey] = val;
	} else {
		it->second = val;
	}
	
	if (isIntra) {
		m_TotalInteractions++;
	}
}

struct interaction_RowSum_worker_t {
	interaction_RowSum_worker_t(int num, int* mps, vector<map<SegmentKey, Interaction > >& m, vector<uint32_t >& rs)
	: nRows(num), midpoints(mps), matrix(m), rowSums(rs)
	{}
	
	int minDelta;
	int maxDelta;
	
	int nRows;
	int* midpoints;
	vector<map<SegmentKey, Interaction > >& matrix;
	vector<uint32_t >& rowSums;
};


// [opt19.2: interaction_RowSum_worker {int}cell-dist,count-LT,UT]
static void interaction_RowSum_worker(void *data, int rowKey, int tid)
{
	interaction_RowSum_worker_t *w = (interaction_RowSum_worker_t*)data;
	
	uint32_t sum=0;
	map<SegmentKey, Interaction > &rowInteractions = w->matrix[rowKey];
	for (map<SegmentKey, Interaction >::iterator it = rowInteractions.begin(); it != rowInteractions.end(); ++it) {
		int delta=w->midpoints[it->first]-w->midpoints[rowKey];
		if(delta>=w->minDelta and delta<=w->maxDelta) {
			sum += it->second.count;
			
			__sync_fetch_and_add(&(w->rowSums[it->first]), it->second.count);
		}
	}
	
	__sync_fetch_and_add(&(w->rowSums[rowKey]), sum);
}

void CCCMatrixInteraction::getRowSums(vector<uint32_t > &res) {
	res.clear();
	res.resize(orderedSegments.size());
	interaction_RowSum_worker_t w((int)orderedSegments.size(), orderedSegmentMidPoints, matrix, res);
	w.maxDelta = maxExpDelta; w.minDelta = minExpDelta;
	kt_for(nThreads, interaction_RowSum_worker, &w, (int)orderedSegments.size());
}


// TODO: optimize
// [opt19.2: getDeltasCount {int,perm}cell-dist,mask-UT]
unsigned long CCCMatrixInteraction::getDeltasCount(int minDelta /* = 0 */, int maxDelta /*=MAXDELTA*/, bool masking /* =false */)
{
	unsigned long size = 0;
	
	for(SegmentKey rowKey=0; rowKey<orderedSegments.size(); ++rowKey) {
		int32_t rowMidPoint = orderedSegmentMidPoints[rowKey];
		map<SegmentKey, Interaction > &rowInteractions = matrix[rowKey];
		for(SegmentKey colKey=rowKey+1; colKey<orderedSegments.size(); ++colKey) {
			map<SegmentKey, Interaction >::iterator it = matrix[rowKey].find(colKey);
			if (it==rowInteractions.end()) {
				int delta=orderedSegmentMidPoints[colKey]-rowMidPoint;
				if(delta>=minDelta and delta<=maxDelta) {
					size++;
				}
			} else {
				if(not (masking and it->second.mask)) {
					int delta=orderedSegmentMidPoints[colKey]-rowMidPoint;
					if(delta>=minDelta and delta<=maxDelta) {
						size++;
					}
				}
			}
		}
	}
	
	return size;
}

//[opt16]
struct getDeltasCount_worker_t {
	getDeltasCount_worker_t(int num, int* mps, vector<map<SegmentKey, Interaction > >& m, DeltaIndicatorCounter* r, unsigned long& n)
	: nRows(num), midpoints(mps), matrix(m), res(r), totalDeltas(n)
	{}
	
	int minDelta;
	int maxDelta;
	bool masking;
	
	int nRows;
	int* midpoints;
	vector<map<SegmentKey, Interaction > >& matrix;
	DeltaIndicatorCounter* res;
	unsigned long& totalDeltas;
};

//[opt16]
// [opt19.2: getDeltasCount_worker {perm}cell-dist-UT]
static void getDeltasCount_worker(void *data, int rowKey, int tid)
{
	getDeltasCount_worker_t *w = (getDeltasCount_worker_t*)data;
	unsigned long totalDeltas = 0;
	
	int32_t rowMidPoint = w->midpoints[rowKey];
	for(SegmentKey colKey=rowKey+1; colKey<w->nRows; ++colKey) {
		int delta=w->midpoints[colKey]-rowMidPoint;
		if(delta>=w->minDelta and delta<=w->maxDelta) {
			__sync_fetch_and_add((w->res+delta), 1);
			totalDeltas++;
		}
	}
	
	__sync_fetch_and_add(&(w->totalDeltas), totalDeltas);
}

// [opt19.2: getMaskedDeltasCount_worker {int,perm}cell-dist,mask-UT]
static void getMaskedDeltasCount_worker(void *data, int rowKey, int tid)
{
	getDeltasCount_worker_t *w = (getDeltasCount_worker_t*)data;
	unsigned long totalDeltas = 0;
	
	int32_t rowMidPoint = w->midpoints[rowKey];
	map<SegmentKey, Interaction > &rowInteractions = w->matrix[rowKey];
	for(SegmentKey colKey=rowKey+1; colKey<w->nRows; ++colKey) {
		int delta=w->midpoints[colKey]-rowMidPoint;
		if(delta>=w->minDelta and delta<=w->maxDelta) {
			map<SegmentKey, Interaction >::iterator it = w->matrix[rowKey].find(colKey);
			if (it==rowInteractions.end()) {
				__sync_fetch_and_add((w->res+delta), 1);
				totalDeltas++;
			} else if (0==it->second.mask) {
				__sync_fetch_and_add((w->res+delta), 1);
				totalDeltas++;
			}
		}
	}
	
	__sync_fetch_and_add(&(w->totalDeltas), totalDeltas);
}

//[opt16]
unsigned long CCCMatrixInteraction::getDeltasCount(DeltaIndicatorCounter* res, int minDelta /* = 0 */, int maxDelta /*=MAXDELTA*/, bool masking /* =false */)
{
	unsigned long totalDeltas = 0;
	getDeltasCount_worker_t w((int)orderedSegments.size(), orderedSegmentMidPoints, matrix, res, totalDeltas);
	w.minDelta = minDelta; w.maxDelta = maxDelta; w.masking = masking;
	if (masking) {
		kt_for(nThreads, getMaskedDeltasCount_worker, &w, (int)orderedSegments.size());
	} else {
		kt_for(nThreads, getDeltasCount_worker, &w, (int)orderedSegments.size());
	}
	
	return totalDeltas;
}

struct getN_worker_t {
	getN_worker_t(int num, int* mps, vector<map<SegmentKey, Interaction > >& m, uint32_t* r)
	: nRows(num), midpoints(mps), matrix(m), res(r)
	{}
	
	int minDelta;
	int maxDelta;
	
	int nRows;
	int* midpoints;
	vector<map<SegmentKey, Interaction > >& matrix;
	uint32_t* res;
};

// [opt19.2: getN_worker {int}cell-dist,count-UT]
static void getN_worker(void *data, int rowKey, int tid)
{
	getN_worker_t *w = (getN_worker_t*)data;
	uint32_t res = 0;

	int32_t rowMidPoint = w->midpoints[rowKey];
	map<SegmentKey, Interaction > &rowInteractions = w->matrix[rowKey];
	for (map<SegmentKey, Interaction >::iterator it = rowInteractions.begin(); it != rowInteractions.end(); ++it) {
		//WCH: NO more mirroring
		// sum a triangle, not the square
		int delta=w->midpoints[it->first]-rowMidPoint;
		if(delta>=w->minDelta and delta<=w->maxDelta) {
			res += it->second.count;
		}
	}
	
	__sync_fetch_and_add(w->res, res);
}

// [opt19.1]
void CCCMatrixInteraction::calculate_cache_N(int minDelta/*=0*/, int maxDelta/*=MAXDELTA*/) {
	if (CACHE_FLAG_N!=(CACHE_FLAG_N&cacheFlag)) {
		uint32_t res = 0;
		getN_worker_t w((int)orderedSegments.size(), orderedSegmentMidPoints, matrix, &res);
		w.minDelta = minDelta; w.maxDelta = maxDelta;
		kt_for(nThreads, getN_worker, &w, (int)orderedSegments.size());
		cachedN = res;
		
		cacheFlag|=CACHE_FLAG_N;
	}
}


struct getStatisticsPerDelta_worker_t {
	getStatisticsPerDelta_worker_t(int num, int* mps, vector<map<SegmentKey, Interaction > >& m, QuantileMapper* pqtm, map<int, DeltaStatistics >& r)
	: nRows(num), midpoints(mps), matrix(m), pQuantileMapper(pqtm), res(r)
	{}
	
	int minDelta;
	int maxDelta;
	
	int nRows;
	int* midpoints;
	vector<map<SegmentKey, Interaction > >& matrix;
	QuantileMapper* pQuantileMapper;
	map<int, DeltaStatistics >& res;
};

// [opt19.2: getN_worker {int}cell-dist,count-UT]
static void getStatisticsPerDeltaNonMasked_worker(void *data, int rowKey, int tid)
{
	getStatisticsPerDelta_worker_t *w = (getStatisticsPerDelta_worker_t*)data;
	
	int32_t rowMidPoint = w->midpoints[rowKey];
	map<SegmentKey, Interaction > &rowInteractions = w->matrix[rowKey];
	for(SegmentKey colKey=rowKey+1; colKey<w->nRows; ++colKey) {
		int delta=w->midpoints[colKey]-rowMidPoint;
		if(delta>=w->minDelta and delta<=w->maxDelta) {
			map<SegmentKey, Interaction >::iterator it = rowInteractions.find(colKey);
			if (it==rowInteractions.end()) {
				__sync_fetch_and_add(&(w->res[w->pQuantileMapper->getDeltaQuantile(delta)].count), 1);
			} else {
				if (!it->second.mask) {
					map<int, DeltaStatistics >::iterator itDeltaStat = w->res.find(w->pQuantileMapper->getDeltaQuantile(delta));
					__sync_fetch_and_add(&(itDeltaStat->second.sum), it->second.count);
					__sync_fetch_and_add(&(itDeltaStat->second.count), 1);
				}
			}
		}
	}
}

static void getStatisticsPerDelta_worker(void *data, int rowKey, int tid)
{
	getStatisticsPerDelta_worker_t *w = (getStatisticsPerDelta_worker_t*)data;
	
	int32_t rowMidPoint = w->midpoints[rowKey];
	map<SegmentKey, Interaction > &rowInteractions = w->matrix[rowKey];
	for(SegmentKey colKey=rowKey+1; colKey<w->nRows; ++colKey) {
		int delta=w->midpoints[colKey]-rowMidPoint;
		if(delta>=w->minDelta and delta<=w->maxDelta) {
			map<SegmentKey, Interaction >::iterator it = rowInteractions.find(colKey);
			if (it==rowInteractions.end()) {
				__sync_fetch_and_add(&(w->res[w->pQuantileMapper->getDeltaQuantile(delta)].count), 1);
			} else {
				map<int, DeltaStatistics >::iterator itDeltaStat = w->res.find(w->pQuantileMapper->getDeltaQuantile(delta));
				__sync_fetch_and_add(&(itDeltaStat->second.sum), it->second.count);
				__sync_fetch_and_add(&(itDeltaStat->second.count), 1);
			}
		}
	}
}

// TODO: OPTIMIZE : multithreading
// [opt19.2: getStatisticsPerDelta {int,perm}cell-dist,mask,count-UT]
void CCCMatrixInteraction::getStatisticsPerDelta(map<int, DeltaStatistics >& res, QuantileMapper* pQuantileMapper, int minDelta/*=8000*/, int maxDelta/*=MAXDELTA*/, bool masking/*=false*/)
{
	getStatisticsPerDelta_worker_t w((int)orderedSegments.size(), orderedSegmentMidPoints, matrix, pQuantileMapper, res);
	w.minDelta = minDelta; w.maxDelta = maxDelta;
	if (masking) {
		kt_for(nThreads, getStatisticsPerDeltaNonMasked_worker, &w, (int)orderedSegments.size());
	} else {
		kt_for(nThreads, getStatisticsPerDelta_worker, &w, (int)orderedSegments.size());
	}
}

struct interaction_Positives_worker_t {
	interaction_Positives_worker_t(int* n, vector<map<SegmentKey, Interaction > >& m)
	: matrix(m), res(n)
	{}
	
	double alpha;
	int cutoff;
	
	vector<map<SegmentKey, Interaction > >& matrix;
	int* res;
};

// [opt19.2: interaction_Positives_A_worker {int}cell-pvalue-UT]
static void interaction_Positives_A_worker(void *data, int rowKey, int tid)
{
	interaction_Positives_worker_t *w = (interaction_Positives_worker_t*)data;
	int count = 0;
	map<SegmentKey, Interaction > &rowInteractions = w->matrix[rowKey];
	for (map<SegmentKey, Interaction >::iterator it = rowInteractions.begin(); it != rowInteractions.end(); ++it) {
		// WCH : NO mirroring, so this check is redundant
		// sum a triangle, not the square
		if (it->second.pvalue <= w->alpha) count++;
	}
	__sync_fetch_and_add(w->res, count);
}

// [opt19.2: interaction_Positives_AC_worker {int}cell-pvalue,count-UT]
static void interaction_Positives_AC_worker(void *data, int rowKey, int tid)
{
	interaction_Positives_worker_t *w = (interaction_Positives_worker_t*)data;
	int count = 0;
	map<SegmentKey, Interaction > &rowInteractions = w->matrix[rowKey];
	for (map<SegmentKey, Interaction >::iterator it = rowInteractions.begin(); it != rowInteractions.end(); ++it) {
		// WCH : NO mirroring, so this check is redundant
		// sum a triangle, not the square
		if (it->second.pvalue <= w->alpha and it->second.count >= w->cutoff)  count++;
	}
	__sync_fetch_and_add(w->res, count);
}

// NO caller
int CCCMatrixInteraction::getNumberOfPositives(double alpha)
{
	assert(isIntra);
	if (0==alpha) return 0;
	
	int res = 0;
	interaction_Positives_worker_t w(&res, matrix);
	w.alpha = alpha;
	kt_for(nThreads, interaction_Positives_A_worker, &w, (int)orderedSegments.size());
	return res;
}

int CCCMatrixInteraction::getNumberOfPositives(double alpha, int cutoff)
{
	assert(isIntra);
	if (0==alpha) return 0;
	
	int res = 0;
	interaction_Positives_worker_t w(&res, matrix);
	w.alpha = alpha; w.cutoff = cutoff;
	kt_for(nThreads, interaction_Positives_AC_worker, &w, (int)orderedSegments.size());
	return res;
}


struct maskByAlpha_worker_t {
	maskByAlpha_worker_t(int num, vector<map<SegmentKey, Interaction > >& m, double a, int c, vector<unsigned long>& tm)
	: nRows(num), matrix(m), alpha(a), cutoff(c), masked(tm)
	{}
	
	int nRows;
	vector<map<SegmentKey, Interaction > >& matrix;
	double alpha;
	int cutoff;
	vector<unsigned long>& masked;
};

static void maskByAlpha_worker(void *data, int rowKey, int tid)
{
	maskByAlpha_worker_t *w = (maskByAlpha_worker_t*)data;
	
	map<SegmentKey, Interaction > &rowInteractions = w->matrix[rowKey];
	for (map<SegmentKey, Interaction >::iterator it = rowInteractions.begin(); it != rowInteractions.end(); ++it) {
		// WCH: NO more mirroring
		// sum a triangle, not the square
		if (it->second.pvalue <= w->alpha and it->second.count >= w->cutoff) {
			it->second.mask = 1;
			w->masked[tid]++;
		}
	}
}

// OPTIMIZE: not a hotspot, skip optimization for now
// [opt19.2: maskByAlpha: {int}cell-pvalue,count,mask-UT]
unsigned long CCCMatrixInteraction::maskByAlpha(double alpha, int cutoff, bool dump/*=false*/)
{
	unsigned long masked = 0;
	vector<unsigned long> tmpRes; tmpRes.resize(nThreads, 0);
	maskByAlpha_worker_t w((int)orderedSegments.size(), matrix, alpha, cutoff, tmpRes);
	kt_for(nThreads, maskByAlpha_worker, &w, (int)orderedSegments.size());
	for(unsigned int i=0; i<nThreads; ++i) { masked += tmpRes[i]; }
	return masked;
}

struct calculate_Pvalues_worker_t {
	calculate_Pvalues_worker_t(int num, int* mps, vector<map<SegmentKey, Interaction > >& m, vector< float>& ge, vector<uint32_t >& gm, spline* sp, vector<double >& rs)
	: nRows(num), midpoints(mps), matrix(m), getExp(ge), getMarginal(gm), delta2ExpSpline(sp), rowSums(rs)
	{}
	
	int n;
	
	int minDelta;
	int maxDelta;
	
	int N2x;
	float L2x;
	
	int nRows;
	int* midpoints;
	vector<map<SegmentKey, Interaction > >& matrix;
	vector< float>& getExp;
	vector<uint32_t >& getMarginal;
	spline* delta2ExpSpline;
	vector<double >& rowSums;
};

// [opt19.2: calculate_Pvales_worker {int}cell-dist,pvalue-UT]
static void calculate_Pvalues_worker(void *data, int rowKey, int tid)
{
	calculate_Pvalues_worker_t *w = (calculate_Pvalues_worker_t*)data;

	// WCH: only need to compute for cell with readings, the rest are defaulted to 1.0
	int rowMidPoint = w->midpoints[rowKey];
	map<SegmentKey, Interaction > &rowInteractions = w->matrix[rowKey];
	for (map<SegmentKey, Interaction >::iterator it = rowInteractions.begin(); it != rowInteractions.end(); ++it) {
		SegmentKey colKey = it->first;
		int delta=w->midpoints[colKey] - rowMidPoint;
		if(delta >= w->minDelta and delta <= w->maxDelta) {
			float Li = w->getExp[rowKey];
			float Lj = w->getExp[colKey];
			int Ni = w->getMarginal[rowKey];
			int Nj = w->getMarginal[colKey];
			
			int Nij = it->second.count;
			float Lij = max((float)w->delta2ExpSpline->cachedAt(delta), fSmallValue); // Arbitrary cutoff, but expectation cannot be negative or zero!

			float omegaDenominator = ((Li-Lij)*(Lj-Lij));
			if (omegaDenominator > fSmallValue ) { // Kind of arbitrary cutoff, but needs to be set at some low value to avoid Omega = Inf
				float Omega = (Lij*(w->L2x-Li-Lj+Lij)) / omegaDenominator;
#ifdef ORIGINAL_CODE
				it->second.pvalue = getPvalNCHG(Nij, Ni,Nj, w->N2x, Omega, fSmallValue, dMinOmega);
#else
				it->second.pvalue = getPvalNCHGNew(Nij, Ni,Nj, w->N2x, Omega);
#endif
			}
			else {
				it->second.pvalue = 1.0;
			}
		} else {
			it->second.pvalue = 1.0;
		}
	}
}

// multi-threading version
void CCCMatrixInteraction::calculatePvalues(int minDelta/*=8000*/, int maxDelta/*=MAXDELTA*/)
{
	assert(isIntra);
	vector<double> rowSums; rowSums.resize(orderedSegments.size());
	calculate_Pvalues_worker_t w((int)orderedSegments.size(), orderedSegmentMidPoints, matrix, m_getExp, m_getMarginal, pDelta2ExpSpline, rowSums);
	w.minDelta = minDelta; w.maxDelta=maxDelta;
	w.N2x = 2 * cachedN; w.L2x = 2.0 * cachedExpectN;
	kt_for(nThreads, calculate_Pvalues_worker, &w, (int)orderedSegments.size());
}

// NO optimization needed
void CCCMatrixInteraction::printPositives(string& chr, double alpha, int cutoff, bool printNij, bool printAll/*=false*/)
{
	for (SegmentKey rowKey=0; rowKey<orderedSegments.size(); ++rowKey) {
		map<SegmentKey, Interaction > &rowInteractions = matrix[rowKey];
		for (map<SegmentKey, Interaction >::iterator it = rowInteractions.begin(); it != rowInteractions.end(); ++it) {
			SegmentKey colKey = it->first;
			//WCH: NO more mirroring
			bool significant = (it->second.pvalue <= alpha and it->second.count >= cutoff);
			if (significant || printAll) {
				cout << chr << "\t" << orderedSegments[rowKey].start << "\t" << orderedSegments[rowKey].end << "\t" <<  chr << "\t" << orderedSegments[colKey].start << "\t" << orderedSegments[colKey].end << "\t" << it->second.pvalue;
				if(printNij) {
					float expect = getExpectElement(rowKey, colKey);
					cout << "\t" << it->second.count << "\t" << ((expect<0) ? 0. : expect);
				}
				
				if (printAll) {
					cout << "\t" << ((significant) ? "*" : ".");
				}
				cout << endl;
			}
		}
	}
}


void CCCMatrixInteraction::setDeltaToExpectationMapper(spline& smoothDelta2ExpMapper, int minDelta/*=8000*/, int maxDelta/*=MAXDELTA*/)
{
	assert(isIntra);
	pDelta2ExpSpline = &smoothDelta2ExpMapper;
	minExpDelta = minDelta;
	maxExpDelta = maxDelta;

	calculate_cache_N(minDelta, maxDelta);

	cacheFlag &= (~CACHE_FLAG_EXPECTN);
	cachedExpectN = 0.0;
	calculate_cache_ExpectN(minDelta, maxDelta);
	
	// cache commonly used values
	cacheRowSums();
	
}

struct getExpectN_worker_t {
	getExpectN_worker_t(int num, int* mps, vector<map<SegmentKey, Interaction > >& m, spline* sp, vector<double>& tr, vector<uint32_t>& tc)
	: nRows(num), midpoints(mps), matrix(m), delta2ExpSpline(sp), tmpRes(tr), tmpCount(tc)
	{}
	
	int minDelta;
	int maxDelta;
	
	int nRows;
	int* midpoints;
	vector<map<SegmentKey, Interaction > >& matrix;
	spline* delta2ExpSpline;
	vector<double>& tmpRes;
	vector<uint32_t>& tmpCount;
};

// [opt19.2: getExpectN_worker {perm}cell-dist-UT]
static void getExpectN_worker(void *data, int rowKey, int tid)
{
	getExpectN_worker_t *w = (getExpectN_worker_t*)data;
	double res = 0.0;
	uint32_t tc = 0;
	
	int32_t rowMidPoint = w->midpoints[rowKey];
	for (SegmentKey colKey=rowKey+1; colKey<w->nRows; ++colKey) {
		int delta=w->midpoints[colKey]-rowMidPoint;
		if(delta>=w->minDelta and delta<=w->maxDelta) {
			res += w->delta2ExpSpline->cachedAt(delta);
			tc++;
		}
	}
	
	w->tmpRes[tid] += res;
	w->tmpCount[tid] += tc;
}

// OPTIMIZE multi-threading
void CCCMatrixInteraction::calculate_cache_ExpectN(int minDelta/*=0*/, int maxDelta/*=MAXDELTA*/)
{
	uint32_t count = 0;
	if (CACHE_FLAG_EXPECTN!=(CACHE_FLAG_EXPECTN&cacheFlag)) {
		if (nThreads>1) {
			vector<double> tmpRes; tmpRes.resize(nThreads, 0.0);
			vector<uint32_t> tmpCount; tmpCount.resize(nThreads, 0);
			getExpectN_worker_t w((int)orderedSegments.size(), orderedSegmentMidPoints, matrix, pDelta2ExpSpline, tmpRes, tmpCount);
			w.minDelta = minDelta; w.maxDelta = maxDelta;
			kt_for(nThreads, getExpectN_worker, &w, (int)orderedSegments.size());
			double dres = 0.0;
			for(unsigned int i=0; i<nThreads; ++i) { dres += tmpRes[i]; count += tmpCount[i]; }
			cachedExpectN = dres;
		}  else {
			float res = 0.0;
			for (SegmentKey rowKey=0; rowKey<orderedSegments.size(); ++rowKey) {
				for (SegmentKey colKey=rowKey+1; colKey<orderedSegments.size(); ++colKey) {
					int delta=orderedSegmentMidPoints[colKey]-orderedSegmentMidPoints[rowKey];
					if(delta>=minDelta and delta<=maxDelta) {
						res += pDelta2ExpSpline->cachedAt(delta);
					}
				}
			}
			cachedExpectN = res;
		}
		
		cacheFlag|=CACHE_FLAG_EXPECTN;
	}
}

void CCCMatrixInteraction::printExpectationMatrix(ostream& os,bool hideMask /*=false*/) {
	// print column header
	os << "Expect";
	for(SegmentKey segmentCounter=0; segmentCounter<orderedSegments.size(); ++segmentCounter) {
		os << " " << segmentCounter;
	}
	os << endl;
	
	// re-written to exclude abs(...) call
	for (SegmentKey rowKey=0; rowKey<orderedSegments.size(); ++rowKey) {
		// print row header
		os << rowKey;
		for (SegmentKey colKey=0; colKey<rowKey; ++colKey) {
			int delta=orderedSegmentMidPoints[rowKey]-orderedSegmentMidPoints[colKey];
			if(delta>=minExpDelta and delta<=maxExpDelta) {
				float res = pDelta2ExpSpline->cachedAt(delta);
				os << " " << res;
			} else {
				os << " 0";
			}
		}
		for (SegmentKey colKey=rowKey; colKey<orderedSegments.size(); ++colKey) {
			int delta=orderedSegmentMidPoints[colKey]-orderedSegmentMidPoints[rowKey];
			if(delta>=minExpDelta and delta<=maxExpDelta) {
				float res = pDelta2ExpSpline->cachedAt(delta);
				os << " " << res;
			} else {
				os << " 0";
			}
		}
		os << endl;
	}
}

struct expectation_RowSum_worker_t {
	expectation_RowSum_worker_t(int num, int* mps, vector<map<SegmentKey, Interaction > >& m, spline* s, vector<float >& rs)
	: nRows(num), midpoints(mps), matrix(m), pDelta2ExpSpline(s), rowSums(rs)
	{}
	
	int n;
	
	int minDelta;
	int maxDelta;
	
	int nRows;
	int* midpoints;
	vector<map<SegmentKey, Interaction > >& matrix;
	spline* pDelta2ExpSpline;
	vector<float >& rowSums;
};

// [opt19.2: expectation_RowSum_worker {perm}cell-dist-LT,UT]
static void expectation_RowSum_worker(void *data, int rowKey, int tid)
{
	expectation_RowSum_worker_t *w = (expectation_RowSum_worker_t*)data;
	
	float sum = 0.0f;
	int rowMidPoint = w->midpoints[rowKey];
	for(SegmentKey colKey = 0; colKey<rowKey; ++colKey) {
		int delta=rowMidPoint - w->midpoints[colKey];
		// WCH: if the full matrix version, these delta constraints will have affected value
		if(delta>=w->minDelta and delta<=w->maxDelta)
			sum += w->pDelta2ExpSpline->cachedAt(delta);
	}
	for(SegmentKey colKey = rowKey+1; colKey<w->nRows; ++colKey) {
		int delta=w->midpoints[colKey]-rowMidPoint;
		// WCH: if the full matrix version, these delta constraints will have affected value
		if(delta>=w->minDelta and delta<=w->maxDelta)
			sum += w->pDelta2ExpSpline->cachedAt(delta);
	}
	w->rowSums[rowKey] = sum;
}

void CCCMatrixInteraction::getExpectRowSums(vector<float > &res)
{
	res.clear();
	res.resize(orderedSegments.size());
	expectation_RowSum_worker_t w((int)orderedSegments.size(), orderedSegmentMidPoints, matrix, pDelta2ExpSpline, res);
	w.maxDelta = maxExpDelta; w.minDelta = minExpDelta;
	kt_for(nThreads, expectation_RowSum_worker, &w, (int)orderedSegments.size());
}


// multithreading version
struct interaction_FalsePositives_worker_t {
	interaction_FalsePositives_worker_t(int num, int* mps, vector<map<SegmentKey, Interaction > >& m, vector< float>& ge, vector<uint32_t >& gm, spline* sp, vector<double >& rs)
	: nRows(num), midpoints(mps), matrix(m), getExp(ge), getMarginal(gm), delta2ExpSpline(sp), rowSums(rs)
	{}
	
	int n;
	
	double alpha;
	int cutoff;
	int minDelta;
	int maxDelta;
	
	int N2x;
	float L2x;

	CCCMatrixInteraction* pMatrix;
	int nRows;
	int* midpoints;
	vector<map<SegmentKey, Interaction > >& matrix;
	vector< float>& getExp;
	vector<uint32_t >& getMarginal;
	spline* delta2ExpSpline;
	vector<double >& rowSums;
};

// [opt19.2: interaction_FalsePositives_worker {perm}cell-dist-UT]
static void interaction_FalsePositives_worker(void *data, int rowKey, int tid)
{
	interaction_FalsePositives_worker_t *w = (interaction_FalsePositives_worker_t*)data;
	
	double res = 0.0;
#ifdef ORIGINAL_CODE
	int cut;
	float Omega = 1.0;
	int quant = -1;
#endif
	int32_t rowMidPoint = w->midpoints[rowKey];
	for (SegmentKey colKey=rowKey+1; colKey<w->nRows; ++colKey) {
		int delta=w->midpoints[colKey]-rowMidPoint;
#ifdef ORIGINAL_CODE
		if(delta >= w->minDelta and delta <= w->maxDelta) {
			float Li = w->getExp[rowKey];
			float Lj = w->getExp[colKey];
			int Ni = w->getMarginal[rowKey];
			int Nj = w->getMarginal[colKey];
			
			float Lij = max((float)w->delta2ExpSpline->cachedAt(delta), fSmallValue);
			
			if (((Li-Lij)*(Lj-Lij)) > 1E-20 ) { // Arbitrary cutoff to avoid Omega = Inf
				Omega = (Lij*(w->L2x-Li-Lj+Lij)) / ((Li-Lij)*(Lj-Lij));
				quant = getQuantileNCHG(1-w->alpha, Ni, Nj, w->N2x, Omega, 1E-50, true);
				cut = max(w->cutoff, 1+quant);
			}
			else {
				cut = 0;
			}
			if(cut > 0 && quant != -1) { // quant == -1 indicates memory issues
				res += getPvalNCHG(cut, Ni, Nj, w->N2x, Omega, dOmegaPrecision, dMinOmega);
			}
		}
#else
		if(delta >= w->minDelta and delta <= w->maxDelta) {
			float Li = w->getExp[rowKey];
			float Lj = w->getExp[colKey];
			int Ni = w->getMarginal[rowKey];
			int Nj = w->getMarginal[colKey];
			
			float Lij = max((float)w->delta2ExpSpline->cachedAt(delta), fSmallValue);
			
			float omegaNumerator = (Lij*(w->L2x-Li-Lj+Lij));
			float omegaDenominator = ((Li-Lij)*(Lj-Lij));
			if (omegaNumerator>=0. && omegaDenominator > fSmallValue ) { // Arbitrary cutoff to avoid Omega = Inf
				float Omega = omegaNumerator / omegaDenominator;
				if (Omega<dMinOmega) {
					// need to perform separately
					int quant = getQuantileNCHG(1-w->alpha, Ni, Nj, w->N2x, Omega);
					if (quant != -1) { // quant == -1 indicates memory issues
						int cut = max(w->cutoff, 1+quant);
						res += getPvalNCHGNew(cut, Ni, Nj, w->N2x, Omega);
					}
				} else {
					// can re-use the FNCHG table
					double pval = getQuantilePvalNCHG(w->cutoff, 1-w->alpha, Ni, Nj, w->N2x, Omega);
					if (pval>0) res += pval;
				}
			}
		}
#endif
	}
	w->rowSums[tid] += res;
}

double CCCMatrixInteraction::estimateFalsePositives(double alpha, int cutoff, int minDelta/*=8000*/, int maxDelta/*=MAXDELTA*/)
{
	assert(isIntra);
	if(alpha == 0) {
		return 0;
	}
	
	vector<double> tmpSums;
	tmpSums.resize(nThreads);
	interaction_FalsePositives_worker_t w((int)orderedSegments.size(), orderedSegmentMidPoints, matrix, m_getExp, m_getMarginal, pDelta2ExpSpline, tmpSums);
	w.alpha = alpha; w.cutoff = cutoff; w.minDelta = minDelta; w.maxDelta=maxDelta;
	w.N2x = 2*cachedN; w.L2x = 2*cachedExpectN;
	kt_for(nThreads, interaction_FalsePositives_worker, &w, (int)orderedSegments.size());
	double res = accumulate(tmpSums.begin(), tmpSums.end(), 0.0);
	
	return res;
}

uint32_t CCCMatrixInteraction::getSegmentCount()
{
	// FIXME: data type size
	return (uint32_t) orderedSegments.size();
}


void CCCMatrixInteraction::cacheRowSums()
{
	getExpectRowSums(m_getExp);
	getRowSums(m_getMarginal);
}

unsigned long CCCMatrixInteraction::getInteractionCount()
{
	return m_TotalInteractions;
}

void CCCMatrixInteraction::setThreads(int n)
{
	nThreads = n;
}

int CCCMatrixInteraction::getWidestInteractionSpan()
{
	if (orderedSegments.size()>0) {
		return orderedSegmentMidPoints[orderedSegments.size()-1] - orderedSegmentMidPoints[0] + 1;
	} else {
		return 0;
	}
}

