CPP = g++
CPPFLAGS = -I . -O3

OBJECTS1 = bait-fisher-helper.o bait-fisher.o Csequence_cluster_and_center_sequence.o mydir-unix.o
OBJECTS2 = bait-filter.o global-types-and-parameters.o

%.o: %.cpp
	$(CPP) $(CPPFLAGS) -c $< -o $@

default: BaitFisher-v1.2.7 BaitFilter-v1.0.5

BaitFisher-v1.2.7: $(OBJECTS1)
	$(CPP) $(CPPFLAGS) -lstdc++ $(OBJECTS1) -o $@

BaitFilter-v1.0.5: $(OBJECTS2)
	$(CPP) $(CPPFLAGS) -lstdc++ $(OBJECTS2) -o $@

clean:
	rm -f *.o
	rm -f BaitFisher-v1.2.7
	rm -f BaitFilter-v1.0.5
