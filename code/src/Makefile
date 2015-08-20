all: experiment 


experiment: run_experiment.cpp gp.cpp kernel_function.cpp prob_util.cpp
	g++ -Wall -std=c++11 $^ -o $@  -I ../inc/

test_gp: test_gp.cpp gp.cpp kernel_function.cpp prob_util.cpp
	g++ -Wall -std=c++11 $^ -o $@  -I ../inc/


clean:
	rm *.o
