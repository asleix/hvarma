all: obj/gradient.so

obj/gradient.so: obj/gradient.o
	gcc -shared -o obj/gradient.so obj/gradient.o

obj/gradient.o:
	gcc -c -O2 -Wall -Werror -fpic src/gradient.c -o obj/gradient.o
	
clean:
	rm obj/gradient.so obj/gradient.o
