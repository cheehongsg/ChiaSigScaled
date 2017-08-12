// ChiaSig: detection of significant interactions in ChIA-PET data using Non-central hypergeometric distribution (NCHG)
// This file is subject to the terms and conditions defined in file 'LICENSE', which is part of this source code package.
// Copyright (2014) Jonas Paulsen


#include "CCCStatistics.h"
#include <iostream>
#include <fstream>

#include <assert.h>

#include "Segment.h"

using namespace std;

extern "C" {
	double cputime();
	double realtime();
}

// TODO: WCH [pending]
void printPositives(string& chr, CCCMatrixInteraction& cpemat, double alpha, int cutoff, bool printNij /*=false*/, bool printAll /*=false*/) {
	assert(cpemat.isIntra);
	cpemat.printPositives(chr, alpha, cutoff, printNij, printAll);
}


// Utils:
// WCH: loop does not terminate
template<typename T> vector<T> getSequence(T from, T to, int length) {
	assert(to >= from);
	vector<T> res;
	T stepSize=(to-from)/(length-1);
	T i = from;
	for(int c=0; c<length; c++) {
		res.push_back(i);
		i+=stepSize;
	}
	return res;
}

pair<double, double> getAlphaRange(map<double, double> alpha2FDR, double minFDR) {
	typedef map<double, double>::iterator it2;
	//WCH
	double alpha=0, FDR=0;
	pair<double, double> res;
	res.first = 0;
	res.second = 0;
	for(it2 iter = alpha2FDR.begin(); iter != alpha2FDR.end(); iter++) {
		alpha = iter->first;
		FDR = iter->second;
		fprintf(stderr, "[M::%s:%d] Re-estimating alpha/FDR=%.6e, FDR=%.6e\n", __func__, __LINE__, alpha, FDR);
		if(FDR > minFDR) {
			if(res.first != 0 ) {
				res.second = alpha;
				fprintf(stderr, "[M::%s:%d] success! FDR(%.6e)>minFDR(%.6e)\n", __func__, __LINE__, FDR, minFDR);
				return res; // success! minFDR is within the range
			}
			else { // minFDR is not within the range. C
				res.first = 0;
				res.second = alpha;
				fprintf(stderr, "[M::%s:%d] minFDR(%.6e) is not within range\n", __func__, __LINE__, minFDR);
				return res; // the right alpha must be somewhere between 0 and minFDR
			}
		}
		res.first = alpha;
	}
	// At this point, the right alpha must be somewhere between alpha and minFDR
	
	fprintf(stderr, "[M::%s:%d] right alpha is between alpha(%.6e) and minFDR(%.6e)\n", __func__, __LINE__, alpha, minFDR);
	res.second = minFDR;
	return res;
}

// TODO: WCH [coded]
static int G_nThreads = 1;
void setThreads(int n) {
	G_nThreads = n;
}

void kt_for(int n_threads, void (*func)(void*,int,int), void *data, int n);

bool chromSegmentCountDesc(const pair<int, uint32_t >& a, const pair<int, uint32_t >& b) {
	return b.second < a.second;
}

// this routine will estimate the largest alpha such that the FDR<fdr
map<double, double>  estimateFDRAlpha(vector<string>& chromosomes, vector<pair<int, uint32_t > >&chromosomesRows, map<ChromosomeIndexType, CCCMatrixInteraction>& cmpmat, double fdr, int refinesteps, int cutoff /*=3*/, int deltaCutoff /*=8000*/, int maxDelta /*=MAXDELTA*/) {
	
	map<double, pair<int, double> > tmp;
	map<double, double> res;
	
	double t_diff;
	double t_start = realtime();

	double alphaLower = 0.0;
	double alphaUpper = fdr / 100.0;
	
	double alpha = alphaUpper;
	int numStep = 1;
	
	// attempt to establish possible range to start iteration
	double t_start_iter = realtime();
	int numPositives = 0;
	double numFalsePositives = 0.0;
	
	for(vector<pair<int, uint32_t > >::iterator itChrom = chromosomesRows.begin(); itChrom!=chromosomesRows.end(); ++itChrom) {
		double t_chrStart = realtime();
		numPositives += cmpmat[itChrom->first].getNumberOfPositives(alpha, cutoff);
		numFalsePositives += cmpmat[itChrom->first].estimateFalsePositives(alpha, cutoff, deltaCutoff, maxDelta);
		t_diff = realtime() - t_chrStart;
		fprintf(stderr, "[M::%s:%d] estimateFDR, chr=%s took %.2f min..\n", __func__, __LINE__, chromosomes[itChrom->first].c_str(), t_diff/60.0);
	}
	
	tmp[alpha] = make_pair(numPositives, numFalsePositives);
	double alphaFDR = (0==numPositives) ? 0.0 : numFalsePositives / numPositives;
	res[alpha] = alphaFDR;
	
	t_diff = realtime() - t_start_iter;
	fprintf(stderr, "[M::%s:%d] step:%d alpha:%.6e FDR:%.6e FP:%.2f P:%u rtmin:%.2f%s\n", __func__, __LINE__, numStep, alpha, alphaFDR, numFalsePositives, numPositives, t_diff/60.0, ((alphaFDR<fdr)?" (*)":""));

	if (alphaFDR >= fdr) {
		// we still have much larger FDR, so we need a more stringent alpha
		alphaUpper = alpha;
		alpha = 0.5 * (alphaLower + alphaUpper);
	} else if (alphaFDR < fdr) {
		// there are more positives that we can report, relax the alpha
		alphaLower = alpha;
		alphaUpper = fdr; // we set the Upper bound to a higher alpha
		alpha = 0.5 * (alphaLower + alphaUpper);
	}
	
	int maxAlphaPositives = 0;
	bool maxAlphaPositivesIncreasing = true;
	numStep++;
	while(numStep<=refinesteps && maxAlphaPositivesIncreasing) {
		
		double t_start_iter = realtime();
		
		int numPositives = 0;
		double numFalsePositives = 0.0;
		
		for(vector<pair<int, uint32_t > >::iterator itChrom = chromosomesRows.begin(); itChrom!=chromosomesRows.end(); ++itChrom) {
			double t_chrStart = realtime();
			numPositives += cmpmat[itChrom->first].getNumberOfPositives(alpha, cutoff);
			numFalsePositives += cmpmat[itChrom->first].estimateFalsePositives(alpha, cutoff, deltaCutoff, maxDelta);
			t_diff = realtime() - t_chrStart;
			fprintf(stderr, "[M::%s:%d] estimateFDR, chr=%s took %.2f min..\n", __func__, __LINE__, chromosomes[itChrom->first].c_str(), t_diff/60.0);
		}
		
		tmp[alpha] = make_pair(numPositives, numFalsePositives);
		double alphaFDR = (0==numPositives) ? 0.0 : numFalsePositives / numPositives;
		res[alpha] = alphaFDR;
		
		t_diff = realtime() - t_start_iter;
		fprintf(stderr, "[M::%s:%d] step:%d alpha:%.6e FDR:%.6e FP:%.2f P:%u rtmin:%.2f%s\n", __func__, __LINE__, numStep, alpha, alphaFDR, numFalsePositives, numPositives, t_diff/60.0, ((alphaFDR<fdr)?" (*)":""));
		
		// let's decide on the next alpha to try
		if (alphaFDR >= fdr) {
			// we still have much larger FDR, so we need a more stringent alpha
			alphaUpper = alpha;
			alpha = 0.5 * (alphaLower + alphaUpper);
		} else if (alphaFDR < fdr) {
			// there are more positives that we can report, relax the alpha
			alphaLower = alpha;
			alpha = 0.5 * (alphaLower + alphaUpper);
			// early termination
			if (numPositives>maxAlphaPositives) {
				maxAlphaPositives = numPositives;
			} else {
				maxAlphaPositivesIncreasing = false;
				if (numStep<refinesteps)
					fprintf(stderr, "[M::%s:%d] expecting insignificant FDR improvement, terminating..\n", __func__, __LINE__);
			}
		}
		
		numStep++;
	}
	
	t_diff = realtime() - t_start;
	fprintf(stderr, "[M::%s:%d] estimateFDR took %.2f min..\n", __func__, __LINE__, t_diff/60.0);
	
	return res;
}


map<double, double>  estimateFDR(vector<string>& chromosomes, vector<pair<int, uint32_t > >&chromosomesRows, map<ChromosomeIndexType, CCCMatrixInteraction>& cmpmat, double alphaFrom, double alphaTo, double alphaN, int cutoff /*=3*/, int deltaCutoff /*=8000*/, int maxDelta /*=MAXDELTA*/) {
	
	//  vector<string> chromosomes=DataReader.getChromosomes();
	map<double, pair<int, double> > tmp;
	map<double, double> res;
	
	assert(alphaFrom <= alphaTo);
	vector<double> alphas = getSequence(alphaFrom, alphaTo, alphaN);
	
	double t_diff;
	double t_start = realtime();
	
	for(vector<pair<int, uint32_t > >::iterator itChrom = chromosomesRows.begin(); itChrom!=chromosomesRows.end(); ++itChrom) {
		string chr=chromosomes[itChrom->first];
		
		double t_chrStart = realtime();
		for(unsigned int j=0; j!=alphas.size(); ++j) {
			double alpha = alphas[j];
			tmp[alpha].first += cmpmat[itChrom->first].getNumberOfPositives(alpha, cutoff);
			tmp[alpha].second += cmpmat[itChrom->first].estimateFalsePositives(alpha, cutoff, deltaCutoff, maxDelta);
		}
		
		t_diff = realtime() - t_chrStart;
		fprintf(stderr, "[M::%s:%d] estimateFDR, chr=%s took %.2f min..\n", __func__, __LINE__, chr.c_str(), t_diff/60.0);
	}
	
	t_diff = realtime() - t_start;
	fprintf(stderr, "[M::%s:%d] estimateFDR took %.2f min..\n", __func__, __LINE__, t_diff/60.0);
	
	for(map<double, pair<int, double> >::iterator iter = tmp.begin(); iter != tmp.end(); iter++) {
		double FDR = (iter->second.first == 0) ? 0 : (iter->second.second) / (iter->second.first);
		res[iter->first] = FDR;
		fprintf(stderr, "[M::%s:%d] alpha:%.6e FDR:%.6e FP:%.2f P:%u \n", __func__, __LINE__, iter->first, FDR, iter->second.second, iter->second.first);
		//cout << "alpha:" << iter->first << " FDR:" << FDR << endl;
	}
	return res;
}


struct estimateDeltaSpline_worker_t {
	estimateDeltaSpline_worker_t(vector<string>& c, map<ChromosomeIndexType, CCCMatrixInteraction>& m, vector<unsigned long>& n, vector<unsigned long>& o, vector<int>& d, vector< map<int, DeltaStatistics > >& ds, QuantileMapper* pq, vector<pair<int, uint32_t > > &ps)
	: chromosomes(c), mat(m), nums(n), offsets(o), deltas(d), deltasStats(ds), pQuantileMapper(pq), parameterSpace(ps)
	{}
	
	int minDelta;
	int maxDelta;
	bool masking;
	
	int n;
	
	vector<string>& chromosomes;
	map<ChromosomeIndexType, CCCMatrixInteraction>& mat;
	vector<int> &deltas;
	vector<pair<int, uint32_t > >& parameterSpace;
	QuantileMapper* pQuantileMapper;
	
	
	vector<unsigned long>& nums;
	vector<unsigned long>& offsets;
	vector<map<int, DeltaStatistics > >& deltasStats;
};

static void estimateDeltaSpline_worker_DeltasCounting(void *data, int i, int tid)
{
	estimateDeltaSpline_worker_t *w = (estimateDeltaSpline_worker_t*)data;
	unsigned int nChrIdx = w->parameterSpace[i].first;
	w->nums[nChrIdx] = w->mat[nChrIdx].getDeltasCount(w->minDelta, w->maxDelta, w->masking);
}

void estimateDeltaSpline(spline& s, vector<double>& sx, vector<string>& chromosomes, vector<pair<int, uint32_t > >&chromosomesRows, map<ChromosomeIndexType, CCCMatrixInteraction>& cmpmat, int minDelta /*=8000*/, int maxDelta /*=MAXDELTA*/, bool masking /*=false*/, int percentiles /*=1000*/) {

	double t_diff;
	// Get all deltas for all chromosomes:
	QuantileMapper quantileMapper;
	// deltas allocations
	int nMaxSpan = 0;
	{
		for(unsigned int i=0; i<chromosomes.size(); ++i) {
			int maxChromSpan = cmpmat[i].getWidestInteractionSpan();
			if (maxChromSpan>nMaxSpan) nMaxSpan = maxChromSpan;
		}
		nMaxSpan++; // 0-based index, so need 1 more allocation
		//cerr << "Largest span " << nMaxSpan << endl;
	}
	// NOTE: deltas stores count, the index is actually the span
	DeltaIndicatorCounter* pDeltas = (DeltaIndicatorCounter*) calloc(nMaxSpan+1, sizeof(DeltaIndicatorCounter));

	double t_start = realtime();
	unsigned long ulTotalDeltas = 0;
	for(unsigned int i=0; i<chromosomes.size(); ++i) {
		double t_subStart = realtime();
		int nChrIdx = chromosomesRows[i].first;
		unsigned long ulDeltas = cmpmat[nChrIdx].getDeltasCount(pDeltas, minDelta, maxDelta, masking);
		ulTotalDeltas += ulDeltas;
		t_diff = realtime() - t_subStart;
		fprintf(stderr, "[M::%s:%d] #deltas(%s)=%lu collected, %.2f min..\n", __func__, __LINE__, chromosomes[nChrIdx].c_str(), ulDeltas, t_diff/60.0);
	}
	t_diff = realtime() - t_start;
	fprintf(stderr, "[M::%s:%d] #deltas(genome)=%lu collected, %.2f min..\n", __func__, __LINE__, ulTotalDeltas, t_diff/60.0);

	double t_subStart = realtime();
	quantileMapper.computeDeltaCountQuantile(ulTotalDeltas, pDeltas, percentiles, masking);
	t_diff = realtime() - t_subStart;
	fprintf(stderr, "[M::%s:%d] delta quantiles computation deltas took %.2f min..\n", __func__, __LINE__, t_diff/60.0);
	free(pDeltas);
	
	t_start = realtime();
	map<int, DeltaStatistics > deltasStats;
	{
		vector<int>& quantileDeltas = quantileMapper.getQuantileDeltas();
		for(vector<int>::iterator it=quantileDeltas.begin(); it!=quantileDeltas.end(); ++it) {
			int deltaQuantile = quantileMapper.getDeltaQuantile(*it);
			deltasStats[deltaQuantile] = DeltaStatistics();
		}
	}
	for(unsigned int i=0; i<chromosomesRows.size(); i++) {
		double t_subStart = realtime();
		int nChrIdx = chromosomesRows[i].first;
		cmpmat[nChrIdx].getStatisticsPerDelta(deltasStats, &quantileMapper, minDelta, maxDelta, masking);
		t_diff = realtime() - t_subStart;
		fprintf(stderr, "[M::%s:%d] delta quantiles statistics computation (%s), %.2f min..\n", __func__, __LINE__, chromosomes[nChrIdx].c_str(), t_diff/60.0);
	}
	t_diff = realtime() - t_start;
	fprintf(stderr, "[M::%s:%d] delta quantiles statistics computation (genome) took %.2f min..\n", __func__, __LINE__, t_diff/60.0);
	
	sx.clear();
	vector<double> exp_y; // y-object for spline generation
	for(map<int, DeltaStatistics>::iterator iter = deltasStats.begin(); iter != deltasStats.end(); iter++) {
		sx.push_back((double)(iter->first));
		exp_y.push_back((double)(iter->second.getMean()));
	}
	
	s.set_points(sx,exp_y);    // currently it is required that X is already sorted
}

void estimateResources(string& filename, vector<string>& chromosomes, vector<pair<int, uint32_t > >&chromosomesRows, map<ChromosomeIndexType, CCCMatrixInteraction>& cmpmat, int minDelta/*=8000*/, int maxDelta/*=MAXDELTA*/, bool masking/*=false*/)
{
	double t_diff;
	double t_start = realtime();
	
	// Get all deltas for all chromosomes:
	vector<int> deltas; // not really used
	vector<map<int, DeltaStatistics > > chromDeltasStats;
	chromDeltasStats.resize(chromosomes.size());
	
	vector<unsigned long> chromDeltaCounts;
	vector<unsigned long> chromDeltaOffsets;
	chromDeltaCounts.resize(chromosomes.size());
	chromDeltaOffsets.resize(chromosomes.size());
	estimateDeltaSpline_worker_t w(chromosomes, cmpmat, chromDeltaCounts, chromDeltaOffsets, deltas, chromDeltasStats, NULL, chromosomesRows);
	w.minDelta = minDelta; w.maxDelta = maxDelta; w.masking = masking;
	kt_for(G_nThreads, estimateDeltaSpline_worker_DeltasCounting, &w, (int)chromosomesRows.size());
	
	unsigned long nTotalDeltas = 0;
	{
		unsigned long nOffset = 0;
		for(map<ChromosomeIndexType, CCCMatrixInteraction>::iterator itChrom=cmpmat.begin(); itChrom!=cmpmat.end(); ++itChrom) {
			nTotalDeltas += chromDeltaCounts[itChrom->first];
			chromDeltaOffsets[itChrom->first] = nOffset;
			nOffset = nTotalDeltas;
		}
	}
	t_diff = realtime() - t_start;
	fprintf(stderr, "[M::%s:%d] Allocating %lu units for deltas, %.2f min\n", __func__, __LINE__, nTotalDeltas, t_diff/60.0);
	
	// loop as table grid!!!
	{
		string ofile(filename); ofile.append(".resource.xls");
		ofstream fileDump;
		fileDump.open(ofile.c_str());
		
		unsigned long totalSegments = 0;
		unsigned long totalInteractions = 0;
		unsigned long totalDeltas = 0;
		fileDump << "#filename\t" << filename.c_str() << endl;
		fileDump << "#chr\tanchors\tinteractions\tdeltas" << endl;
		for(unsigned int i=0; i<chromosomes.size(); ++i) {
			totalSegments += cmpmat[i].getSegmentCount();
			totalInteractions += cmpmat[i].getInteractionCount();
			totalDeltas += chromDeltaCounts[i];
			fileDump << chromosomes[i];
			fileDump << "\t" << cmpmat[i].getSegmentCount();
			fileDump << "\t" << cmpmat[i].getInteractionCount();
			fileDump << "\t" << chromDeltaCounts[i];
			fileDump << endl;
		}
		fileDump << "GRAND TOTAL";
		fileDump << "\t" << totalSegments;
		fileDump << "\t" << totalInteractions;
		fileDump << "\t" << totalDeltas;
		fileDump << endl;
		
		fileDump.close();
	}
}

void getChromosomesDataSizeDesc(vector<string>& chromosomes, map<ChromosomeIndexType, CCCMatrixInteraction>& cmpmat, vector<pair<int, uint32_t > >& chromosomesRows)
{
	chromosomesRows.clear();
	chromosomesRows.reserve(chromosomes.size());
	for(unsigned int i=0; i<chromosomes.size(); ++i) {
		chromosomesRows.push_back(make_pair(i, cmpmat[i].getSegmentCount()));
	}
	sort(chromosomesRows.begin(), chromosomesRows.end(), chromSegmentCountDesc);
}

void sweepQunatilesDeltaSpline(string& filename, vector<string>& chromosomes, vector<pair<int, uint32_t > >&chromosomesRows, map<ChromosomeIndexType, CCCMatrixInteraction>& cmpmat, vector<unsigned int >& percentilesList, int minDelta/*=8000*/, int maxDelta/*=MAXDELTA*/, bool masking/*=false*/)
{
	double t_diff;
	
	// Get all deltas for all chromosomes:
	vector<map<int, DeltaStatistics > > chromDeltasStats;
	chromDeltasStats.resize(chromosomes.size());
	
	std::sort(percentilesList.begin(), percentilesList.end(), std::greater<unsigned int>());
	
	// deltas allocations
	int nMaxSpan = 0;
	{
		for(unsigned int i=0; i<chromosomes.size(); ++i) {
			int maxChromSpan = cmpmat[i].getWidestInteractionSpan();
			if (maxChromSpan>nMaxSpan) nMaxSpan = maxChromSpan;
		}
	}
	DeltaIndicatorCounter* pDeltas = (DeltaIndicatorCounter*) calloc(nMaxSpan+1, sizeof(DeltaIndicatorCounter));
	
	vector<vector<splineType > > cells;
	cells.resize(percentilesList.size());
	
	int dmax = (maxDelta == MAXDELTA) ? nMaxSpan / 100 : maxDelta;
	vector<double> d = getSequence((double)minDelta, (double)dmax, 2000);
	for(unsigned int percentileIdx = 0; percentileIdx<percentilesList.size(); ++percentileIdx) {
		cells[percentileIdx].resize(d.size());
	}
	
	
	double t_start = realtime();
	unsigned long ulTotalDeltas = 0;
	for(unsigned int i=0; i<chromosomes.size(); ++i) {
		double t_subStart = realtime();
		int nChrIdx = chromosomesRows[i].first;
		unsigned long ulDeltas = cmpmat[nChrIdx].getDeltasCount(pDeltas, minDelta, maxDelta, masking);
		ulTotalDeltas += ulDeltas;
		t_diff = realtime() - t_subStart;
		fprintf(stderr, "[M::%s:%d] #deltas(%s)=%lu collected, %.2f min..\n", __func__, __LINE__, chromosomes[nChrIdx].c_str(), ulDeltas, t_diff/60.0);
	}
	t_diff = realtime() - t_start;
	fprintf(stderr, "[M::%s:%d] #deltas(genome)=%lu collected, %.2f min..\n", __func__, __LINE__, ulTotalDeltas, t_diff/60.0);

	for(unsigned int percentileIdx = 0; percentileIdx<percentilesList.size(); ++percentileIdx) {
		double t_subStart = realtime();
		QuantileMapper qtm;
		spline s;
		qtm.computeDeltaCountQuantile(ulTotalDeltas, pDeltas, percentilesList[percentileIdx], masking);
		t_diff = realtime() - t_subStart;
		fprintf(stderr, "[M::%s:%d] delta quantiles(%u) computation deltas took %.2f min..\n", __func__, __LINE__, percentilesList[percentileIdx], t_diff/60.0);
		
		t_subStart = realtime();
		map<int, DeltaStatistics > deltasStats;
		{
			vector<int>& quantileDeltas = qtm.getQuantileDeltas();
			for(vector<int>::iterator it=quantileDeltas.begin(); it!=quantileDeltas.end(); ++it) {
				int deltaQuantile = qtm.getDeltaQuantile(*it);
				deltasStats[deltaQuantile] = DeltaStatistics();
			}
		}
		for(unsigned int i=0; i<chromosomesRows.size(); i++) {
			double t_subStart = realtime();
			int nChrIdx = chromosomesRows[i].first;
			cmpmat[nChrIdx].getStatisticsPerDelta(deltasStats, &qtm, minDelta, maxDelta, masking);
			t_diff = realtime() - t_subStart;
			fprintf(stderr, "[M::%s:%d] delta quantiles statistics computation (%s), %.2f min..\n", __func__, __LINE__, chromosomes[nChrIdx].c_str(), t_diff/60.0);
		}
		t_diff = realtime() - t_subStart;
		fprintf(stderr, "[M::%s:%d] delta quantiles(%u) statistics computation deltas took %.2f min..\n", __func__, __LINE__, percentilesList[percentileIdx], t_diff/60.0);
		
		vector<double> delta_x; // x-object for spline generation
		vector<double> exp_y; // y-object for spline generation
		for(map<int, DeltaStatistics>::iterator iter = deltasStats.begin(); iter != deltasStats.end(); iter++) {
			delta_x.push_back((double)(iter->first));
			exp_y.push_back((double)(iter->second.getMean()));
		}
		
		// If x is not sorted, an error will occur.
		s.set_points(delta_x,exp_y);    // currently it is required that X is already sorted

		// keep results
		for(unsigned int i=0; i < d.size(); i++) {
			cells[percentileIdx][i] = s.cachedAt(d[i]) ;
		}
	}
	
	free(pDeltas);
	
	// write result in a table grid
	if (percentilesList.size()>0)
	{
		string ofile(filename); ofile.append(".quantiles.sweep.xls");
		ofstream fileDump;
		fileDump.open(ofile.c_str());
		
		// print header
		fileDump << "delta";
		for(unsigned int percentileIdx = 0; percentileIdx<percentilesList.size(); ++percentileIdx) {
			fileDump << "\tsmoothed-" << percentilesList[percentileIdx] ;
		}
		fileDump << endl;

		for(unsigned int i=0; i < d.size(); i++) {
			fileDump << d[i];
			for(unsigned int percentileIdx = 0; percentileIdx<percentilesList.size(); ++percentileIdx) {
				fileDump << "\t" << cells[percentileIdx][i];
			}
			fileDump << endl;
		}
		fileDump.close();
	}
}


void computeInteractionsPvalues(vector<string>& chromosomes, vector<pair<int, uint32_t > >&chromosomesRows, map<ChromosomeIndexType, CCCMatrixInteraction>& cmpmat, spline& s, int deltaCutoff, int minDelta/*=8000*/, int maxDelta/*=MAXDELTA*/)
{
	// [exp{chromsome}, pval{chromosome}] is faster than [chromosome{exp, pval}]
	for(unsigned int j=0; j<chromosomesRows.size(); j++) {
		int nChrIdx = chromosomesRows[j].first;
		string chr=chromosomes[nChrIdx];
		double t_loopStart = realtime();
		cmpmat[nChrIdx].setDeltaToExpectationMapper(s, minDelta, MAXDELTA);
		double t_diff = realtime() - t_loopStart;
		fprintf(stderr, "[M::%s:%d] %s, createExpectationMatrix took %.2f min..\n", __func__, __LINE__, chr.c_str(), t_diff/60.0);
	}
	for(unsigned int j=0; j<chromosomesRows.size(); j++) {
		int nChrIdx = chromosomesRows[j].first;
		string chr=chromosomes[nChrIdx];
		double t_loopStart = realtime();
		cmpmat[nChrIdx].calculatePvalues(deltaCutoff, MAXDELTA);// During refinement, all P-values for all deltas are considered!
		
		double t_diff = realtime() - t_loopStart;
		fprintf(stderr, "[M::%s:%d] %s, calculatePvalues took %.2f min..\n", __func__, __LINE__, chr.c_str(), t_diff/60.0);
	}
}


// Add the possible template classes explicitly, to allow linker to find them:
// These are the only allowed:
template vector<double> getSequence(double from, double to, int length);
template vector<int> getSequence(int from, int to, int length);
template vector<float> getSequence(float from, float to, int length);
