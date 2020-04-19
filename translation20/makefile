tred: main.o buildk.o tracek.o oneiteration.o
	g++ -std=c++11 -g -o tred main.o buildk.o tracek.o oneiteration.o 

main.o: main.cpp buildk.h oneiteration.h parameters.h tracek.h
	g++ -std=c++11 -c -g main.cpp 

buildk.o: buildk.cpp buildk.h parameters.h
	g++ -c -g buildk.cpp

tracek.o: tracek.cpp tracek.h buildk.h parameters.h
	g++ -c -g tracek.cpp

oneiteration.o: oneiteration.cpp oneiteration.h buildk.h tracek.h parameters.h
	g++ -c -g oneiteration.cpp
