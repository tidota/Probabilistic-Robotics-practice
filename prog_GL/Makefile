CC=g++
CFLAGS=-Wall -Ofast -std=c++11
PRG=GLsim
OBJ=MV3.o env_xyz.o
LIB=-lgsl

all: $(PRG)

$(PRG): simumain.cpp $(OBJ) 
	$(CC) $(CFLAGS) -o $@ $^ $(LIB)

%.o: %.cpp %.hpp 
	$(CC) $(CFLAGS) -o $@ -c $< $(LIB) 

clean:
	rm -rf $(OBJ) $(PRG)
