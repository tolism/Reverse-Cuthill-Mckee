CC=gcc
CFLAGS= -O3 -fopenmp  -Wall -std=c99 
SDIR=./src
LDIR=./lib
IDIR=./inc
INCLUDES = -I$(IDIR)
LIBRARIES = -L$(LDIR)
LIBS= $(LDIR)/rcm_sequential.a $(LDIR)/rcm_parallel.a 
EXE= $(SDIR)/rcm_sequential $(SDIR)/rcm_parallel
LDFLAGS = -lm



all: $(LIBS) $(EXE)
lib: $(LIBS)



$(SDIR)/rcm_sequential.o: $(SDIR)/sequential.c
	$(CC) $(CFLAGS) $(LDFLAGS) $(INCLUDES) -o $@ -c $^ 


$(SDIR)/rcm_parallel.o: $(SDIR)/parallel.c
	$(CC) $(CFLAGS) $(LDFLAGS) $(INCLUDES) -o $@ -c $^


$(SDIR)/%.o: $(SDIR)/%.c
	$(CC) $(CFLAGS) $(INCLUDES) -o $@ -c $^

$(LDIR)/%.a: $(SDIR)/%.o $(SDIR)/queue.o
	ar rcs $@ $^


$(LDIR)/rcm_sequential.a: $(SDIR)/rcm_sequential.o $(SDIR)/queue.o
	ar rcs $@ $^ 


$(LDIR)/rcm_parallel.a: $(SDIR)/rcm_parallel.o $(SDIR)/queue.o
	ar rcs $@ $^ 




$(SDIR)/rcm_sequential: $(SDIR)/sequential.c $(LDIR)/rcm_sequential.a
	$(CC) $(CFLAGS) $(INCLUDES) $(LIBRARIES) -o $@ $^ $(LDFLAGS)


$(SDIR)/rcm_parallel: $(SDIR)/parallel.c $(LDIR)/rcm_parallel.a
	$(CC) $(CFLAGS) $(INCLUDES) $(LIBRARIES) -o $@ $^ $(LDFLAGS)





clean:
	$(RM) $(SDIR)/*.o *~ $(LDIR)/rcm_*.a $(EXE)
