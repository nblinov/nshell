CC=g++
CFLAGS= -O3 #-fopenmp #-floop-parallelize-all -ftree-parallelize-loops=4 ##-fopenmp 
DEFINES = -DUSE_OPEN_GL
## for Mac OS X
#GLFLAGS = -framework OpenGL -framework GLUT -lobjc
## For Linux
LDFLAGS = -L/usr/X11R6/lib -lm #-pg
GLFLAGS = -lglut -lGLU -lGL -lXmu -lXi -lXext -lX11 -lm 

SRC=src


radial_collapse: radial_collapse.cpp $(SRC)/shell.hpp $(SRC)/nbody_system.hpp $(SRC)/integrator.hpp
	$(CC) -o radial_collapse radial_collapse.cpp $(CFLAGS) $(LDFLAGS) 

shell_collapse_emd: shell_collapse_emd.cpp $(SRC)/shell.hpp $(SRC)/nbody_system.hpp $(SRC)/nbody_system_emd.hpp $(SRC)/integrator.hpp
	$(CC) -o shell_collapse_emd shell_collapse_emd.cpp $(CFLAGS) $(LDFLAGS)

self_similar_collapse: self_similar_collapse.cpp $(SRC)/shell.hpp $(SRC)/nbody_system.hpp $(SRC)/integrator.hpp
	$(CC) -o self_similar_collapse self_similar_collapse.cpp $(CFLAGS) $(LDFLAGS) 

two_shell: two_shell.cpp $(SRC)/shell.hpp $(SRC)/nbody_system.hpp $(SRC)/integrator.hpp
	$(CC) -o two_shell two_shell.cpp $(CFLAGS) $(LDFLAGS) 

clean:
	rm -f radial_collapse
	rm -f self_similar_collapse
	rm -f two_shell
	rm -f shell_collapse_emd
############################################################
# OLD
solar_openGL: solar_openGL.cpp particle.hpp vec2.hpp
	$(CC) -o solar_openGL solar_openGL.cpp $(CFLAGS) $(DEFINES) $(LDFLAGS) $(GLFLAGS)

shell_clump: shell_clump.cpp shell.hpp
	$(CC) -o shell_clump shell_clump.cpp $(CFLAGS) $(DEFINES) $(LDFLAGS) $(GLFLAGS)

spherical_collapse_open_GL: spherical_collapse.cpp shell.hpp
	$(CC) -o spherical_collapse spherical_collapse.cpp $(CFLAGS) $(DEFINES) $(LDFLAGS)  $(GLFLAGS)

spherical_collapse: spherical_collapse.cpp shell.hpp
	$(CC) -o spherical_collapse spherical_collapse.cpp $(CFLAGS) $(LDFLAGS)  $(GLFLAGS)

one_shell: one_shell.cpp shell.hpp
	$(CC) -o one_shell one_shell.cpp $(CFLAGS) $(LDFLAGS)  $(GLFLAGS)

