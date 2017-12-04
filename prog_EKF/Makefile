CC=g++
CFLAGS=-Wall -std=c++11
PRG=EKFsim
OBJ=MV3.o env_xyz.o

all: $(PRG)

$(PRG): simumain.cpp $(OBJ) 
	$(CC) $(CFLAGS) -o $@ $^

%.o: %.cpp %.hpp 
	$(CC) $(CFLAGS) -o $@ -c $< 

clean:
	rm -rf $(OBJ) $(PRG)
