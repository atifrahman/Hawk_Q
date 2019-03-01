all:
	g++ hawk_q.cpp cpp/src/*.cpp -lm -lpthread -o hawk_q
	g++ hawk.cpp cpp/src/*.cpp -lm -lpthread -o hawk
	g++ bonf_fasta.cpp -o bonf_fasta
	g++ kmersearch.cpp -o kmersearch
	g++ kmersummary.cpp -o kmersummary
	g++ preProcess.cpp -o preProcess
	g++ convertToFasta.cpp -o convertToFasta