INCDIR = /usr/X11/include
XLIBS  = /usr/X11R6/lib

pendulum:		pendulum.mk pendulum.c get_color.ro.c
		gcc -I $(INCDIR) -o pendulum pendulum.c get_color.ro.c \
			-L $(XLIBS) \
			$(EXTLIBS) -lXaw -lXmu -lXt -lX11 -lm
