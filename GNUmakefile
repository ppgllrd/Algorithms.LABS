
CC      	= gcc 
#LINK        = gcc -o$(MAIN)  $(OBJECTS)  -g -lintl -lpthread libs/dcmt0.6.1/lib/libdcmt.a
LINK        = gcc -o$(MAIN)  $(OBJECTS)  -g -rdynamic -lpthread libs/dcmt0.6.1/lib/libdcmt.a
STRIP       = strip


RM      = rm

#CFLAGS          = -g  -gstabs -std=gnu99 -O3 -lpthread -DNDEBUG
CFLAGS          = -g -rdynamic -gstabs -std=gnu99 -O3 -DNDEBUG -D__RUN_IN_LINUX__


OBJECTS     = dynamicMem.o random.o GA.o timer.o LABS.o TabuSearch.o


MAIN        = ma

default     :: $(MAIN)


.c.o:
	$(CC) -c $(CFLAGS) $<


$(MAIN) :	$(OBJECTS)
		$(LINK) 

clean::
		$(RM) *.o
		$(RM) $(MAIN)

