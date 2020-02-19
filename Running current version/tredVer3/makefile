tred: main.o buildk.o tracek.o errorsarray.o printcompact.o createarray.o oneiteration.o
	g++ -g -o tred main.o buildk.o tracek.o errorsarray.o printcompact.o createarray.o oneiteration.o 

main.o: main.cpp errorsarray.h createarray.h oneiteration.h parameters.h
	g++ -c -g main.cpp 

buildk.o: buildk.cpp  buildk.h 
	g++ -c -g buildk.cpp

tracek.o: tracek.cpp tracek.h buildk.h
	g++ -c -g tracek.cpp

errorsarray.o: errorsarray.cpp errorsarray.h
	g++ -c -g errorsarray.cpp 

printcompact.o: printcompact.cpp printcompact.h
	g++ -c -g printcompact.cpp

createarray.o: createarray.cpp createarray.h
	g++ -c -g createarray.cpp

oneiteration.o: oneiteration.cpp oneiteration.h buildk.h tracek.h errorsarray.h printcompact.h parameters.h
	g++ -c -g oneiteration.cpp
