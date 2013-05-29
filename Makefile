
CXX      = g++   
LIBGAB   = /home/gabriel_renaud/lib/

CXXFLAGS = -Wall -lm -O3 -lz -I${LIBGAB}  -c
LDFLAGS  = -lz


all: arrayDesign 

arrayDesign.o:	arrayDesign.cpp
	${CXX} ${CXXFLAGS} arrayDesign.cpp


arrayDesign:	arrayDesign.o ${LIBGAB}utils.o  
	${CXX} -o $@ $^ $(LDLIBS) $(LDFLAGS)

clean :
	rm -f arrayDesign.o arrayDesign

