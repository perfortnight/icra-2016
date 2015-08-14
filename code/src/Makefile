all: test_gp

test_gp: test_gp.cpp gp.cpp kernel_function.cpp
	g++ -Wall -std=c++11 $^ -o $@  -I ../inc/


clean:
	rm *.o
