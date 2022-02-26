CC=g++
CFLAGS=-c
LDFLAGS=-lpthread
SOURCES=main.cpp extra_functions.cpp functions.cpp
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=a

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(OBJECTS) $(LDFLAGS) -o $@

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@

clean:
	rm -f $(OBJECTS)
