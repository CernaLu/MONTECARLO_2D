CC = gcc-6
CFLAGS = -lm -O3

../bin/bidimensional_MC: random.o functions.o f_atmpt.o 
	$(CC) -o ../bin/bidimensional_MC random.o functions.o\
	      f_atmpt.o  $(CFLAGS)

random.o: random.c
	$(CC) -c random.c  $(CFLAGS)

functions.o: functions.c
	$(CC) -c functions.c  $(CFLAGS)

f_atmpt.o: f_atmpt.c
	$(CC) -c f_atmpt.c  $(CFLAGS)

clean:
	rm *.o
