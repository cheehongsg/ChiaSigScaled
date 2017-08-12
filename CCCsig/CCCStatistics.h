// ChiaSig: detection of significant interactions in ChIA-PET data using Non-central hypergeometric distribution (NCHG)
// This file is subject to the terms and conditions defined in file 'LICENSE', which is part of this source code package.
// Copyright (2014) Jonas Paulsen

#ifndef GUARD_CCCStatistics
#define GUARD_CCCStatistics

#include "CCCMatrix.h"
#include "Segment.h"

#include <map>
#include <vector>

#include <algorithm>
#include <numeric> 

#include "../spline/spline.h" // Spline-smoothing

void setThreads(int n);

void printPositives(std::string&, CCCMatrixInteraction&, double, int, bool printNij=false, bool printAll=false);

// Utils:
template<typename T> std::vector<T> getSequence(T, T, int);

std::pair<double, double> getAlphaRange(std::map<double, double> , double);

std::map<double, double>  estimateFDRAlpha(std::vector<std::string>&, std::vector<std::pair<int, uint32_t > >&chromosomesRows, std::map<ChromosomeIndexType, CCCMatrixInteraction>&, double, int, int cutoff=3, int deltaCutoff=8000, int maxDelta=MAXDELTA);

std::map<double, double>  estimateFDR(std::vector<std::string>&, std::vector<std::pair<int, uint32_t > >&chromosomesRows, std::map<ChromosomeIndexType, CCCMatrixInteraction>&, double, double, double, int cutoff=3, int deltaCutoff=8000, int maxDelta=MAXDELTA);

void computeInteractionsPvalues(std::vector<std::string>& chromosomes, std::vector<std::pair<int, uint32_t > >&, std::map<ChromosomeIndexType, CCCMatrixInteraction>& cmpmat, spline& smoothDelta2ExpMapper, int deltaCutoff, int minDelta=8000, int maxDelta=MAXDELTA);

void estimateDeltaSpline(spline&, std::vector<double>&, std::vector<std::string>&, std::vector<std::pair<int, uint32_t > >&, std::map<ChromosomeIndexType, CCCMatrixInteraction>&, int minDelta=8000, int maxDelta=MAXDELTA, bool masking=false, int percentiles=1000);

void estimateResources(std::string&, std::vector<std::string>&, std::vector<std::pair<int, uint32_t > >&, std::map<ChromosomeIndexType, CCCMatrixInteraction>&, int minDelta=8000, int maxDelta=MAXDELTA, bool masking=false);

void getChromosomesDataSizeDesc(std::vector<std::string>& chromosomes, std::map<ChromosomeIndexType, CCCMatrixInteraction>& cmpmat, std::vector<std::pair<int, uint32_t > >& chromosomesRows);

void sweepQunatilesDeltaSpline(std::string&, std::vector<std::string>&, std::vector<std::pair<int, uint32_t > >&, std::map<ChromosomeIndexType, CCCMatrixInteraction>&, std::vector<unsigned int >&, int minDelta=8000, int maxDelta=MAXDELTA, bool masking=false);

#endif
