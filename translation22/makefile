tred: main.o aux_main.o buildk.o tracek.o oneiteration.o NewMatrixL.o
	g++ -std=c++11 -g -o tred main.o aux_main.o buildk.o tracek.o oneiteration.o NewMatrixL.o -pthread -lpthread

main.o: main.cpp aux_main.h mingw.thread.h buildk.h oneiteration.h parameters.h tracek.h
	g++ -std=c++11 -c -g main.cpp

aux_main.o: aux_main.cpp aux_main.h buildk.h oneiteration.h parameters.h tracek.h
	g++ -std=c++11 -c -g aux_main.cpp

buildk.o: buildk.cpp buildk.h parameters.h
	g++ -c -g buildk.cpp

tracek.o: tracek.cpp tracek.h buildk.h parameters.h NewMatrixL.h
	g++ -c -g tracek.cpp

oneiteration.o: oneiteration.cpp oneiteration.h buildk.h tracek.h parameters.h NewMatrixL.h
	g++ -c -g oneiteration.cpp

NewMatrixL.o: NewMatrixL.cpp NewMatrixL.h parameters.h
	g++ -std=c++11 -c -g NewMatrixL.cpp

clean:
	rm -f *.o