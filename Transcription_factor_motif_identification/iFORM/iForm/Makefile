iForm-iForm.o: iForm.c
	cc -c -o iForm-iForm.o iForm.c
	gcc -Wall -Wno-unused -DUNIX -D__USE_FIXED_PROTOTYPES__ -O3 -o iForm iForm-iForm.o  ./lib/libcommon.a -lm  -lgsl -lgslcblas
	rm -r iForm-iForm.o