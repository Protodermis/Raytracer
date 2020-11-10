CFLAGS = -std=c99
OBJS = raytracer.o vecmat.o msg.o
PROG = ray

all: $(PROG)

$(PROG): $(OBJS)
	$(CC) $(CFLAGS) $(OBJS) -o $@ -lm

raytracer.o: vecmat.h msg.h
vecmat.o: vecmat.h msg.h
msg.o:msg.h

clean:
	rm ray *.o *.png
