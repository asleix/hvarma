all: gradient.so

gradient.so: gradient.o
	gcc -shared -o gradient.so gradient.o

gradient.o:
	gcc -c -O2 -Wall -Werror -fpic gradient.c -o gradient.o
	
clean:
	rm gradient.so gradient.o
