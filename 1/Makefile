CC=g++
CFLAGS=-c -O3
LDFLAGS=
SOURCES=main.cpp extra_functions.cpp functions.cpp
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=a 

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@

clean:
	rm $(OBJECTS)
