CFLAGS = -O3 

TCLAPINCLUDE=tclap-1.2.1/include
CCCSIGPATH=CCCsig
SPLINEPATH=spline
STOCCPATH=stocc
UNITTESTPATH=unittest

CHIASIG=$(CCCSIGPATH)/kthread.cpp $(CCCSIGPATH)/CCCDataReader.cpp $(CCCSIGPATH)/CCCMatrix.cpp $(CCCSIGPATH)/Segment.cpp $(CCCSIGPATH)/CCCStatistics.cpp $(CCCSIGPATH)/ChiaSig.cpp $(SPLINEPATH)/spline.cpp $(STOCCPATH)/fnchyppr.cpp $(STOCCPATH)/userintf.cpp 

CHIASIGTARGET=ChiaSig

.PHONY: ChiaSig unittest

ChiaSig:
	g++ $(CFLAGS) -I $(TCLAPINCLUDE) -lpthread $(CHIASIG) -o $(CHIASIGTARGET)

unittest:
	g++ $(CCCSIGPATH)/CCCMatrix.cpp $(CCCSIGPATH)/Segment.cpp $(UNITTESTPATH)/UnitTestMatrix.cpp -o $(UNITTESTPATH)/test_matrix
	g++ $(CCCSIGPATH)/CCCDataReader.cpp $(CCCSIGPATH)/CCCMatrix.cpp $(CCCSIGPATH)/Segment.cpp $(UNITTESTPATH)/UnitTestReadData.cpp -o $(UNITTESTPATH)/test_reader
	g++ $(CCCSIGPATH)/CCCMatrix.cpp $(CCCSIGPATH)/Segment.cpp $(CCCSIGPATH)/CCCStatistics.cpp $(UNITTESTPATH)/UnitTestStats.cpp $(SPLINEPATH)/spline.cpp stocc/fnchyppr.cpp stocc/mersenne.cpp stocc/stoc1.cpp stocc/stoc3.cpp stocc/userintf.cpp stocc/wnchyppr.cpp -o $(UNITTESTPATH)/test_stats
	g++ $(CCCSIGPATH)/CCCDataReader.cpp $(CCCSIGPATH)/CCCMatrix.cpp $(CCCSIGPATH)/Segment.cpp $(CCCSIGPATH)/CCCStatistics.cpp $(UNITTESTPATH)/UnitTestDeltas.cpp $(SPLINEPATH)/spline.cpp $(STOCCPATH)/fnchyppr.cpp $(STOCCPATH)/mersenne.cpp $(STOCCPATH)/stoc1.cpp $(STOCCPATH)/stoc3.cpp $(STOCCPATH)/userintf.cpp $(STOCCPATH)/wnchyppr.cpp -o $(UNITTESTPATH)/test_delta
	g++ $(CCCSIGPATH)/CCCMatrix.cpp $(CCCSIGPATH)/Segment.cpp $(CCCSIGPATH)/CCCStatistics.cpp $(UNITTESTPATH)/UnitTestExpectation.cpp $(SPLINEPATH)/spline.cpp $(STOCCPATH)/fnchyppr.cpp $(STOCCPATH)/mersenne.cpp $(STOCCPATH)/stoc1.cpp $(STOCCPATH)/stoc3.cpp $(STOCCPATH)/userintf.cpp $(STOCCPATH)/wnchyppr.cpp -o $(UNITTESTPATH)/test_expecation
	g++ $(CCCSIGPATH)/CCCDataReader.cpp $(CCCSIGPATH)/CCCMatrix.cpp $(CCCSIGPATH)/Segment.cpp $(CCCSIGPATH)/CCCStatistics.cpp $(UNITTESTPATH)/UnitTestQuantileFunc.cpp $(SPLINEPATH)/spline.cpp $(STOCCPATH)/fnchyppr.cpp $(STOCCPATH)/mersenne.cpp $(STOCCPATH)/stoc1.cpp $(STOCCPATH)/stoc3.cpp $(STOCCPATH)/userintf.cpp $(STOCCPATH)/wnchyppr.cpp -o $(UNITTESTPATH)/test_qfunc
	g++ $(CCCSIGPATH)/CCCDataReader.cpp $(CCCSIGPATH)/CCCMatrix.cpp $(CCCSIGPATH)/Segment.cpp $(CCCSIGPATH)/CCCStatistics.cpp $(UNITTESTPATH)/UnitTestHypTest.cpp $(SPLINEPATH)/spline.cpp $(STOCCPATH)/fnchyppr.cpp $(STOCCPATH)/mersenne.cpp $(STOCCPATH)/stoc1.cpp $(STOCCPATH)/stoc3.cpp $(STOCCPATH)/userintf.cpp $(STOCCPATH)/wnchyppr.cpp -o $(UNITTESTPATH)/test_hyptest
