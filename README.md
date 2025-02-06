<header>

# A double-pendulum simulation

This is a C application written using the X11 libraries that demonstrates
the behavior of a double pendulum, which consists of two rigid rods joined
at their ends, and pivoting around a central point.

The simulation uses equations derived from a Lagrangian model of the
physics of the system to compute the positions and velocities of each
pendulum for the next time step.

The user interface allows the physics to be either real or fake, where the
fake physics simply treats each pendulum separately.

</header>

## Supported environments

The code was written on macOS and requires the Quartz package to be installed
in order to compile and run it.

Go to https://www.xquartz.org to download and install XQuartz for macOS.

The code should also build and run on most versions of Linux that have the
X11 libraries installed.

## Building and running the code

Simply type:

<pre>
make -f pendulum.mk
</pre>

which should create the executable, which can then be run:

<pre>
./pendulum
</pre>

