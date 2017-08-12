// ChiaSig: detection of significant interactions in ChIA-PET data using Non-central hypergeometric distribution (NCHG)
// This file is subject to the terms and conditions defined in file 'LICENSE', which is part of this source code package.
// Copyright (2014) Jonas Paulsen

#define PROGRAM_NAME "ChiaSig"
#define PROGRAM_VERSION "v1.19.44"

#include "../CCCsig/CCCDataReader.h"
#include "../CCCsig/Segment.h"
#include "../CCCsig/CCCMatrix.h"
#include "../CCCsig/CCCStatistics.h"

#include <set>
#include <iostream>
#include <fstream>
#include <sstream>
#include <assert.h>

#include <cstdio>
#include <cstdlib>
#include <vector>

#include "../tclap-1.2.1/include/tclap/CmdLine.h" // Command-line parsing

using namespace std;

// WCH:
#include <sys/time.h>
#include <sys/resource.h>


double G_t_real;

#ifdef __cplusplus
extern "C" {
#endif
	double cputime();
	double realtime();
#ifdef __cplusplus
}
#endif


/*********
 * Timer *
 *********/

double cputime()
{
	struct rusage r;
	getrusage(RUSAGE_SELF, &r);
	return r.ru_utime.tv_sec + r.ru_stime.tv_sec + 1e-6 * (r.ru_utime.tv_usec + r.ru_stime.tv_usec);
}

double realtime()
{
	struct timeval tp;
	struct timezone tzp;
	gettimeofday(&tp, &tzp);
	return tp.tv_sec + tp.tv_usec * 1e-6;
}
// END - WCH:

#define DEFAULT_NQUANT 500

void print_spline(spline s1, spline s2, vector<double> deltas, unsigned int nQuant=DEFAULT_NQUANT, bool errstream = true) {
	(errstream ? cerr : cout) << "#delta" << "\t" << "smoothed-" << nQuant << "\t" << "refined-" << nQuant << endl;
	for(unsigned int i=0; i < deltas.size(); i++) {
		(errstream ? cerr : cout)  << deltas[i] << "\t" << s1.cachedAt(deltas[i]) << "\t" << s2.cachedAt(deltas[i]) << endl;
	}
}

int main(int argc, char** argv) {
	G_t_real = realtime();
	double t_diff;
	
	try {
		TCLAP::CmdLine cmd("ChiaSig is a program to find significant interactions in ChIA-PET datasets, using the Non-central Hypergeometric (NCHG) distribution. For details about this statistical method, see Paulsen et al. 'A statistical model of ChIA-PET data for accurate detection of chromatin 3D interactions', Nucleic acids research. (2014).\n Copyright (2014) Jonas Paulsen\nThis version is the time-and-space optimized version written by Chee-Hong WONG <cheehongsg@gmail.com> based on Paulsen work (v0.93).\n Copyright (2015-6) Chee-Hong WONG" , ' ', "0.93");
		
		TCLAP::SwitchArg resourceEstimationArg("z","resource","Only report essential resource metrices.", false);
		cmd.add( resourceEstimationArg);
		
		TCLAP::MultiArg<unsigned int> quantilesSweepArg("N","Quantiles","Specify (multiple) number of quantiles to use for estimating the expected number of interactions given the genomic distance. (See -n option.) This output a tabulated list of results for comparison to select the best quantiles.", false, "int");
		cmd.add( quantilesSweepArg);
		
		TCLAP::ValueArg<unsigned int> percentilesArg("n","quantiles","Number of quantiles to use for estimating the expected number of interactions given the genomic distance (default is 500). Higher numbers (1000 and more) are appropriate for very large datasets, with many interacting anchors, while lower numbers (some hundred) will be appropriate for smaller datasets.",false,DEFAULT_NQUANT,"int");
		cmd.add( percentilesArg );
		
		TCLAP::SwitchArg skipFDRArg("f","skipFDR","Skip FDR correction, report raw P-values only",false);
		cmd.add( skipFDRArg );
		
		
		TCLAP::ValueArg<double> alphaArg("a","alpha","False discovery rate for selecting significant interactions (default is 0.05).",false,0.05,"double");
		cmd.add( alphaArg );
		
		
		TCLAP::ValueArg<unsigned int> cutoffArg("c","cutoff","Minimum allowed number of interactions for a given pair of anchors, to consider it significant (default is 3)",false,3,"int");
		cmd.add( cutoffArg );
		
		TCLAP::ValueArg<double> refinealphaArg("A","refinealpha","False discovery rate used for refinement (default is 0.01).",false,0.01,"double");
		cmd.add( refinealphaArg );
		
		
		TCLAP::ValueArg<unsigned int> refinecutoffArg("C","refinecutoff","Minimum allowed number of interactions for a given pair of anchors, to consider it significant during refinement (default is 3)",false,3,"int");
		cmd.add( refinecutoffArg );
		
		TCLAP::ValueArg<unsigned int> stepsArg("s","steps","Number of steps used to calculate false discovery rate (default is 22). Higher numbers will give more accurate estimation of the false discovery rate, but will also be slower.",false,22,"int");
		cmd.add( stepsArg );
		
		TCLAP::ValueArg<unsigned int> refinestepsArg("S","refinesteps","Number of steps used to calculate false discovery rate, during refinement (default is 22). Higher numbers will give more accurate estimation of the false discovery rate, but will also be slower.",false,22,"int");
		cmd.add( refinestepsArg );
		
		TCLAP::SwitchArg printdeltaSwitch("d","printdelta","Print estimated expectation values to stderr", false);
		cmd.add( printdeltaSwitch );
		
		
		TCLAP::SwitchArg printInteractionCountsSwitch("p","printcounts","Print observed and expected number of interactions for each significant interactions, in addition to P-values", false);
		cmd.add( printInteractionCountsSwitch );
		TCLAP::SwitchArg printAllInteractionsSwitch("P","printall","Print all interactions regardless of significance", false);
		cmd.add( printAllInteractionsSwitch );
		
		
		TCLAP::SwitchArg onlyprintdeltaSwitch("o","onlyprintdelta","Estimate the expectation values one time (no refinement), then print estimated expectation values to stdout and exit. No P-values are calculated. Suitable for exploring the effect of choosing different number of quantiles (-n).", false);
		cmd.add( onlyprintdeltaSwitch );
		
		
		
		TCLAP::ValueArg<int> mindeltaArg("m","mindelta","Minumum genomic distance (in bp) allowed between anchor pairs, below which interactions are excluded. (default is 8000 bp)",false,8000,"int");
		cmd.add( mindeltaArg );
		
		
		TCLAP::ValueArg<int> maxdeltaArg("M","maxdelta","Maximum genomic distance (in bp) allowed between anchor pairs, above which interactions are excluded. (Default is no maximum value)",false,MAXDELTA,"int");
		cmd.add( maxdeltaArg );
		
		TCLAP::ValueArg<string> selectChrArg("u","onlyChr","Use only the specified chromosome for P-value and FDR calculations. Default behaviour is to use all chromosomes. Even if a chromosome is specified using this argument, all chromosomes will be used for estimation of expected number of interactions.",false,"","string");
		cmd.add( selectChrArg );
		
		
		TCLAP::ValueArg<int> threadArg("t","thread","Number of threads (default is 1)",false,1,"int");
		cmd.add( threadArg );
		
		
		TCLAP::UnlabeledValueArg<string> nolabel( "filename", "Input file of the format chrA startA endA chrB startB endB I(A,B) where I(A,B) gives the number of interactions between A and B. This corresponds to the BEDPE format described here: http://bedtools.readthedocs.org/en/latest/content/general-usage.html#bedpe-format", true, "/dev/null", "filename"  );
		cmd.add( nolabel );
		
		cmd.parse( argc, argv );
		
		// WCH : to facilitate comparison, record the command invocation
		cout << "# version: " << PROGRAM_NAME << " " << PROGRAM_VERSION << endl;
		cout << "# commands:"; for (int i=0; i<argc; ++i) cout << " " << argv[i]; cout << endl;
		
		bool resourceEstimation = resourceEstimationArg.getValue();
		vector<unsigned int> quantilesSweep = quantilesSweepArg.getValue();
		unsigned int percentiles = percentilesArg.getValue();
		int minDelta = mindeltaArg.getValue();
		int maxDelta = maxdeltaArg.getValue();
		bool printDelta = printdeltaSwitch.getValue();
		bool onlyprintDelta = onlyprintdeltaSwitch.getValue();
		bool printAllInteractions = printAllInteractionsSwitch.getValue();
		string fileName = nolabel.getValue();
		double alphaCut = alphaArg.getValue();
		unsigned int cutoff = cutoffArg.getValue();
		int nThreads = threadArg.getValue();
		
		unsigned int steps= stepsArg.getValue();
		unsigned int refinesteps= refinestepsArg.getValue();
		
		double refinealphaCut = refinealphaArg.getValue();
		unsigned int refinecutoff = refinecutoffArg.getValue();
		
		bool skipFDR = skipFDRArg.getValue();
		bool printNij = printInteractionCountsSwitch.getValue();
		string selectChr = selectChrArg.getValue();
		
		assert(cutoff >= 0);
		assert(refinecutoff >= 0);
		assert(steps >= 0);
		assert(refinesteps >= 0);
		assert(minDelta >= 0);
		assert(refinesteps > 0);
		assert(alphaCut > 0 and alphaCut <= 1);
		assert(refinealphaCut > 0 and refinealphaCut <= 1);
		
		setThreads(nThreads);
		
		// READING THE DATA:
		double t_start = realtime();
		int deltaCutoff=minDelta;
		map<ChromosomeIndexType, CCCMatrixInteraction> contactMatrices;
		vector<string> chromosomes;
		vector<pair<int, uint32_t > > sizeOrderedChromosomes;
		{
			CCCDataReader dr(fileName);
			dr.buildContactMatrices();
			dr.getChromosomes(chromosomes);
			t_diff = realtime() - t_start;
			fprintf(stderr, "[M::%s:%d] dr.buildContactMatrices() took %.2f min..\n", __func__, __LINE__, t_diff/60.0);
			
			t_start = realtime();
			// Build contact-matrices:
			assert(0==selectChr.size() || dr.isChromosomePresent(selectChr));
			for(ChromosomeIndexType i=0; i<chromosomes.size(); i++) {
				string chr=chromosomes[i];
				contactMatrices[i] = dr.getContactMatrix(i);
				contactMatrices[i].InteractionLoaded(minDelta, maxDelta);
				contactMatrices[i].setThreads(nThreads);
			}
			t_diff = realtime() - t_start;
			fprintf(stderr, "[M::%s:%d] Build contact-matrices took %.2f min..\n", __func__, __LINE__, t_diff/60.0);
			
			// keep the chromosomes in descending order of data size
			getChromosomesDataSizeDesc(chromosomes, contactMatrices, sizeOrderedChromosomes);
		}
		
		if (resourceEstimation) {
			t_start = realtime();
			fprintf(stderr, "[M::%s:%d] Performing resource estimation..\n", __func__, __LINE__);
			estimateResources(fileName, chromosomes, sizeOrderedChromosomes, contactMatrices, deltaCutoff, MAXDELTA, false);
			t_diff = realtime() - t_start;
			fprintf(stderr, "[M::%s:%d] Resource estimation took %.2f min..\n", __func__, __LINE__, t_diff/60.0);

			t_diff = realtime() - G_t_real;
			fprintf(stderr, "[M::%s:%d] Total time %.2f min..\n", __func__, __LINE__, t_diff/60.0);
			return 0;
		}
		
		if (quantilesSweep.size()>0) {
			t_start = realtime();
			fprintf(stderr, "[M::%s:%d] Performing %lu quantiles sweep..\n", __func__, __LINE__, quantilesSweep.size());
			sweepQunatilesDeltaSpline(fileName, chromosomes, sizeOrderedChromosomes, contactMatrices, quantilesSweep, deltaCutoff, MAXDELTA, false);
			t_diff = realtime() - t_start;
			fprintf(stderr, "[M::%s:%d] %lu quantiles sweep took %.2f min..\n", __func__, __LINE__, quantilesSweep.size(), t_diff/60.0);
			
			t_diff = realtime() - G_t_real;
			fprintf(stderr, "[M::%s:%d] Total time %.2f min..\n", __func__, __LINE__, t_diff/60.0);
			return 0;
		}
		
		// Estimate the spline-object, for the first time:
		t_start = realtime();
		spline s;
		vector<double> deltas;
		estimateDeltaSpline(s, deltas, chromosomes, sizeOrderedChromosomes, contactMatrices, deltaCutoff, MAXDELTA, false, percentiles); // Spline is estimated on all delta>minDelta
		t_diff = realtime() - t_start;
		fprintf(stderr, "[M::%s:%d] estimateDeltaSpline took %.2f min..\n", __func__, __LINE__, t_diff/60.0);
		t_start = realtime();
		t_diff = realtime() - t_start;
		fprintf(stderr, "[M::%s:%d] spline object access took %.2f min..\n", __func__, __LINE__, t_diff/60.0);
		
		if(onlyprintDelta) {
			int dmax = (maxDelta == MAXDELTA) ? deltas.back() / 100 : maxDelta;
			vector<double> d = getSequence((double)deltaCutoff, (double)dmax, 2000);
			print_spline(s, s, d, percentiles, false);
			
			t_diff = realtime() - G_t_real;
			fprintf(stderr, "[M::%s:%d] Total time %.2f min..\n", __func__, __LINE__, t_diff/60.0);
			return 0;
		}

  
		t_start = realtime();
		// Calculating P-values (first time):
		computeInteractionsPvalues(chromosomes, sizeOrderedChromosomes, contactMatrices, s, deltaCutoff, minDelta, MAXDELTA);
		t_diff = realtime() - t_start;
		fprintf(stderr, "[M::%s:%d] p-values calculation first time took %.2f min..\n", __func__, __LINE__, t_diff/60.0);
		
		// IF ONLY RAW P-VALUES SHOULD BE PRINTED::
		if(skipFDR) {
			for(ChromosomeIndexType i=0; i<chromosomes.size(); i++) {
				string chr=chromosomes[i];
				printPositives(chr, contactMatrices[i], 1.0, 0,printNij,printAllInteractions); // Everything is printed. Including 0s.
			}
			
			t_diff = realtime() - G_t_real;
			fprintf(stderr, "[M::%s:%d] Total time %.2f min..\n", __func__, __LINE__, t_diff/60.0);
			return 0;
		}
		
		
		// ----
		
		if (selectChr != "") {
			vector<string> selChr;
			selChr.push_back(selectChr);
			chromosomes = selChr;
		}
		
		
		// Masking, based on P-values, using the double procedure:
		t_start = realtime();
		// In original code, each "steps" will iterate 6 alpha values or 5 intervals
		map<double, double> FDRs = estimateFDRAlpha(chromosomes, sizeOrderedChromosomes, contactMatrices, refinealphaCut, refinesteps, refinecutoff, deltaCutoff, maxDelta);
		t_diff = realtime() - t_start;
		fprintf(stderr, "[M::%s:%d] Masking based on p-values, estimateFDRAlpha took %.2f min..\n", __func__, __LINE__, t_diff/60.0);
		t_start = realtime();
		pair<double, double> alphaRange = getAlphaRange(FDRs, refinealphaCut);
		t_diff = realtime() - t_start;
		fprintf(stderr, "[M::%s:%d] getAlphaRange(alpha:%.6e, FDR:%.6e) took %.2f min..\n", __func__, __LINE__, alphaRange.first, FDRs[alphaRange.first], t_diff/60.0);
		
		cerr << "# Masking P-value cutoff: " << scientific << setprecision(6) << alphaRange.first << ". FDR<=" << scientific << setprecision(6) << FDRs[alphaRange.first] << endl;
		cout << "# Masking P-value cutoff: " << scientific << setprecision(6) << alphaRange.first << ". FDR<=" << scientific << setprecision(6) << FDRs[alphaRange.first] << endl;
		
		// Masking data at the selected range:
		t_start = realtime();
		unsigned long genomeMasked = 0;
		for(ChromosomeIndexType i=0; i<chromosomes.size(); i++) {
			double t_loopStart = realtime();
			string chr=chromosomes[i];
			unsigned long masked = contactMatrices[i].maskByAlpha(alphaRange.first, cutoff, true);
			genomeMasked+=masked;
			t_diff = realtime() - t_loopStart;
			fprintf(stderr, "[M::%s:%d] #masked(%s)=%lu masked in maskByAlpha took %.2f min..\n", __func__, __LINE__, chr.c_str(), masked, t_diff/60.0);
		}
		t_diff = realtime() - t_start;
		fprintf(stderr, "[M::%s:%d]#masked(genome)=%lu took %.2f min..\n", __func__, __LINE__, genomeMasked, t_diff/60.0);
		
		t_start = realtime();
		//fdr:
		// Re-estimating the delta-function:
		spline sRefined;
		vector<double> deltasRefined;
		estimateDeltaSpline(sRefined, deltasRefined, chromosomes, sizeOrderedChromosomes, contactMatrices, deltaCutoff, MAXDELTA, true, percentiles);
		t_diff = realtime() - t_start;
		fprintf(stderr, "[M::%s:%d] estimateDeltaSpline took %.2f min..\n", __func__, __LINE__, t_diff/60.0);
		t_diff = realtime() - t_start;
		fprintf(stderr, "[M::%s:%d] estimateDeltaSpline spline vector took %.2f min..\n", __func__, __LINE__, t_diff/60.0);

		
		// Re-estimating expectations and P-values, now with updated expectations/lambdas:
		t_start = realtime();
		computeInteractionsPvalues(chromosomes, sizeOrderedChromosomes, contactMatrices, sRefined, deltaCutoff, minDelta, MAXDELTA);
		t_diff = realtime() - t_start;
		fprintf(stderr, "[M::%s:%d] Re-estimating expectations and P-values, now with updated expectations/lambdas, loop took %.2f min..\n", __func__, __LINE__, t_diff/60.0);
		
		// Re-estimating FDR on new P-values:
		t_start = realtime();
		// In original code, each "steps" will iterate 6 alpha values or 5 intervals
		FDRs = estimateFDRAlpha(chromosomes, sizeOrderedChromosomes, contactMatrices, alphaCut, steps, cutoff, deltaCutoff, maxDelta);
		t_diff = realtime() - t_start;
		fprintf(stderr, "[M::%s:%d] Re-estimating FDR on new P-values, estimateFDRAlpha took %.2f min..\n", __func__, __LINE__, t_diff/60.0);
		t_start = realtime();
		alphaRange = getAlphaRange(FDRs, alphaCut);
		t_diff = realtime() - t_start;
		fprintf(stderr, "[M::%s:%d] Re-estimating FDR on new P-values, getAlphaRange(alpha:%.6e, FDR:%.6e) took %.2f min..\n", __func__, __LINE__, alphaRange.first, FDRs[alphaRange.first], t_diff/60.0);
		
		cerr << "# Refined P-value cutoff: " << scientific << setprecision(6) << alphaRange.first << ". FDR<=" << scientific << setprecision(6) << FDRs[alphaRange.first] << endl;
		cout << "# Refined P-value cutoff: " << scientific << setprecision(6) << alphaRange.first << ". FDR<=" << scientific << setprecision(6) <<
		FDRs[alphaRange.first] << endl;
		
		//Print the final, significant interactions:
		for(ChromosomeIndexType i=0; i<chromosomes.size(); i++) {
			string chr=chromosomes[i];
			printPositives(chr, contactMatrices[i], alphaRange.first, cutoff, printNij,printAllInteractions);
		}
		
		// Print the estimated expectation values:
		if(printDelta) {
			print_spline(s, sRefined, deltas, percentiles);
		}
		
		
	} catch (TCLAP::ArgException &e)  // catch any exceptions
	{
		std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl;
	}
	
	t_diff = realtime() - G_t_real;
	fprintf(stderr, "[M::%s:%d] Total time %.2f min..\n", __func__, __LINE__, t_diff/60.0);
	
	return 0;
}

