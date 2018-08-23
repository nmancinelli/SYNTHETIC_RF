# Makefile for PropMatrix program

#F77 = f77
F77 = gfortran
FFLAGS = -g
SACLIBRARY = /usr/local/sac/lib/libsac.a
BINDIR = ../bin

synth:		synth.o synth.sub.o
		$(F77)  -o $(BINDIR)/synth synth.o synth.sub.o

synth.o:	synth.f
		$(F77) -g -c -w synth.f

synth.sub.o:	synth.sub.f
		$(F77) -g -c -w synth.sub.f

sourc1:		sourc1.o sourc1.o fork.fold.o
		$(F77)  -o $(BINDIR)/sourc1 sourc1.o fork.fold.o

sourc1.o:	sourc1.f
		$(F77) -std=legacy -g -c -w sourc1.f

fork.fold.o:	fork.fold.f
		$(F77) -g -c -w fork.fold.f

saccpt:	saccpt.f	
		$(F77) $(FFLAGS) -o $(BINDIR)/saccpt saccpt.f $(SACLIBRARY) -ldl

clean:
		trash $(BINDIR)/sourc1
		trash $(BINDIR)/synth
		trash *.o

