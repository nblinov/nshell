CC=g++
CFLAGS= -O3 #-fopenmp #-floop-parallelize-all -ftree-parallelize-loops=4 ##-fopenmp 
DEFINES = -DUSE_OPEN_GL
## for Mac OS X
#GLFLAGS = -framework OpenGL -framework GLUT -lobjc
## For Linux
LDFLAGS = -L/usr/X11R6/lib -lm
GLFLAGS = -lglut -lGLU -lGL -lXmu -lXi -lXext -lX11 -lm 

all: solar_openGL shell_clump spherical_collapse

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

radial_collapse: radial_collapse.cpp shell.hpp nbody_system.hpp integrator.hpp
	$(CC) -o radial_collapse radial_collapse.cpp $(CFLAGS) $(LDFLAGS) 

self_similar_collapse: self_similar_collapse.cpp shell.hpp nbody_system.hpp integrator.hpp
	$(CC) -o self_similar_collapse self_similar_collapse.cpp $(CFLAGS) $(LDFLAGS) 

two_shell: two_shell.cpp shell.hpp nbody_system.hpp integrator.hpp
	$(CC) -o two_shell two_shell.cpp $(CFLAGS) $(LDFLAGS) 
