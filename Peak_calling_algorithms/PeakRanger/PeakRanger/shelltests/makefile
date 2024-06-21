TESTS = $(wildcard ./*.sh)
OBJ=$(TESTS:%.sh=%.sho)


%.sho: %.sh 
	./$< 

all: $(OBJ)



