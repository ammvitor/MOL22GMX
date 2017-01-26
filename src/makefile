#!/bin/bash
CC=g++
CFLAGS= -I.
DEPS=mol2.h gaffparser.h
SOURCES= main.cpp mol2.cpp gaffparser.cpp
OBJECTS= main.o mol2.o gaffparser.o
EXECUTABLE=../MOL22GMX

%.o: %.cpp $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(OBJECTS) -o $@ $(LDFLAGS)

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@

clean:
	rm $(OBJECTS)
