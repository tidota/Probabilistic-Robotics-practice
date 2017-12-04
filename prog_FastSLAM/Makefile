CC=g++
CFLAGS=-Wall -std=c++11 -O2
PRG=FastSLAM1.0
OBJ=MV2.o env_xy.o MVn_gsl.o
LIB=-lgsl -lm -lgslcblas

all: $(PRG)

$(PRG): simumain.cpp $(OBJ) 
	$(CC) $(CFLAGS) -o $@ $^ $(LIB)

%.o: %.cpp %.hpp 
	$(CC) $(CFLAGS) -o $@ -c $< 

clean:
	rm -rf $(OBJ) $(PRG)
