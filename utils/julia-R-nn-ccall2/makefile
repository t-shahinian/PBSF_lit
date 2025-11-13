
all:
	clang++ -c -fpic nn.cpp -o nn.o
	clang++ -shared -o libnn.dylib nn.o
	mkdir -p ../sim/pics

clean:
	rm -f *.o libnn.dylib

