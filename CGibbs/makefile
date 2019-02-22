OUTPUT = CGibbs

CPP = g++

HEADER = energy.h fitness_calculate.h random.h data_table_io.h z_score_map_log10.h parameters.h

INCLUDE = -I.. -I/data/env  # for tnt_126 access

#DEFINE= -DTNT_BOUNDS_CHECK -D_GLIBCXX_DEBUG -g -Wall
#CFLAG = -O0 $(DEFINE)

CFLAG = -O2 # -Wall
#LIB = -L/.../js5084/.../usr/lib

main: $(HEADER) main.o energy.o fitness_calculate.o random.o
	$(CPP) -o $(OUTPUT) $(CFLAG) $(INCLUDE) main.o energy.o fitness_calculate.o random.o 

main.o: $(HEADER) main.cpp
	$(CPP) $(INCLUDE) $(CFLAG) -c main.cpp

energy.o: $(HEADER) energy.cpp
	$(CPP) $(INCLUDE) $(CFLAG) -c energy.cpp

fitness_calculate.o: $(HEADER) fitness_calculate.cpp
	$(CPP) $(INCLUDE) $(CFLAG) -c fitness_calculate.cpp

random.o: $(HEADER) random.cpp
	$(CPP) $(INCLUDE) $(CFLAG) -c random.cpp

clean:
	rm *.o $(OUTPUT)
