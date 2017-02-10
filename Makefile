CXXFLAGS = -g -O2 -std=c++11 -Wall -pipe

gms: gms.o geo.o util.o
	g++ -pthread -std=c++11 $^ -o $@ /usr/local/lib/libs2.a -lcrypto

clean:
	rm -f *.o gms
