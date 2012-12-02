CC=g++
LIBS=-lboost_thread -Iann_1.1.2/include -Lann_1.1.2/lib -lANN
CFLAGS=-c -Wall -O3
LDFLAGS=-Iann_1.1.2/include -Lann_1.1.2/lib -lANN -lboost_thread 
SOURCES=avd.cpp globals.cpp Miniball_dynamic_d.cpp scale.cpp fst.cpp wspd.cpp overlap.cpp representatives.cpp query.cpp 
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=avd

all: $(EXECUTABLE)


$(EXECUTABLE): ann $(OBJECTS)
	$(CC) $(OBJECTS) $(LDFLAGS) -o $@

ann:
	#cd ann_1.1.2 ; make linux-g++
	
.cpp.o:
	$(CC) $(CFLAGS) $(LIBS) $< -o $@

clean:
	rm *.o avd
	#cd ann_1.1.2 ; make clean
