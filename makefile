OBJS = diagnostics.o initial_conditions.o evol.o particle2.o
CC = g++
DEBUG = -g
CFLAGS = -O3 -std=gnu++11  -c -fopenmp $(DEBUG)
LFLAGS = -O3 -std=gnu++11  -fopenmp $(DEBUG)
#CFLAGS = -O3 -std=gnu++11  -c  $(DEBUG)
#LFLAGS = -O3 -std=gnu++11   $(DEBUG)

p1 : $(OBJS)
	$(CC) $(LFLAGS) $(OBJS) -o p1

particle2.o : particle2.cpp initial_conditions.h evol.h diagnostics.h constants.h
	$(CC) $(CFLAGS) particle2.cpp

diagnostics.o : diagnostics.h diagnostics.cpp evol.h constants.h
	$(CC) $(CFLAGS) diagnostics.cpp

initial_conditions.o : initial_conditions.h initial_conditions.cpp constants.h
	$(CC) $(CFLAGS) initial_conditions.cpp

evol.o : evol.h evol.cpp constants.h
	$(CC) $(CFLAGS) evol.cpp

clean:
	\rm *.o p1
