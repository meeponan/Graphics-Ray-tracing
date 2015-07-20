CC = gcc
NAME = liu_xiuwen
LIBS = -lglut32 -lglu32 -lopengl32
CFLAGS = -O3

LIBS =  -L/usr/X11R6/lib/ -O2 -lglut -lGLU -lGL -lXmu -lXt -lSM -lICE -lXext -lX11 -lXi -lXext -lX11 -lm

lab4: lab4.o SSD_util.o
	$(CC) -o lab4 lab4.o SSD_util.o $(LIBS)
lab4_extra: lab1_extra.o SSD_util.o
	$(CC) -o lab4_extra lab4_extra.o SSD_util.o $(LIBS)
.c.o: 
	$(CC)  $(CFLAGS) -c  $(COPT) $<
tar:
	tar cvfz lab4_$(NAME).tar.gz *.c *.h 
	ls -l lab4_$(NAME).tar.gz
run:
	./lab4 lab4_scene1.ssd
	./lab4 lab4_scene2.ssd
	./lab4 lab4_scene3.ssd
	
clean:
	rm  *.o

