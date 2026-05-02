# Coupled Pendulum Simulation

A real-time physics simulation of two spring-coupled pendulums written in C. The system is numerically integrated using the **Velocity Verlet** method and rendered with **raylib**.

---

## Demo

Two pendulums swing from fixed pivots and are connected by a spring at a configurable attachment point. The spring connecting rod is color-coded from **blue → red** based on the instantaneous potential energy stored in the spring, making energy transfer between the two pendulums easy to visualise.

---

## Physics

The system models two simple pendulums of masses $m_1$, $m_2$ and lengths $l_1$, $l_2$ coupled by a spring of constant $k$. The equations of motion are:

$$\ddot{\theta}_1 = -\frac{g}{l_1}\sin\theta_1 - \frac{k \cdot s_1}{m_1 l_1}(s_1\theta_1 - s_2\theta_2)$$

$$\ddot{\theta}_2 = -\frac{g}{l_2}\sin\theta_2 + \frac{k \cdot s_2}{m_2 l_2}(s_1\theta_1 - s_2\theta_2)$$

where $s_1$ and $s_2$ are the spring contact points (fraction of rod length from pivot).

Integration is performed with the **Velocity Verlet** (Störmer–Verlet) scheme, which is symplectic and therefore conserves energy well over long simulations.

---

## Project Structure

```
couple_pendulum/
├── main.c                  # Entry point, rendering loop, user input
├── pendulum/
│   ├── pendulum.h          # Pendulum struct and acceleration computation
│   └── pendulum.c
├── integrator/
│   ├── verlet.h            # Velocity Verlet integrator (half-step scheme)
│   └── verlet.c
└── linearalgebra/
    ├── linearalgebra.h     # Minimal matrix library (create, free)
    └── linearalgebra.c
```

---

## Dependencies

| Library | Purpose |
|---------|---------|
| [raylib](https://www.raylib.com/) | Window creation & 2D rendering |
| `math.h` | `sin`, `cos`, `pow` |

Install raylib on Debian/Ubuntu:
```bash
sudo apt install libraylib-dev
```

Or follow the [raylib build instructions](https://github.com/raysan5/raylib#build-and-installation) for your platform.

---

## Building

```bash
gcc main.c pendulum/pendulum.c integrator/verlet.c linearalgebra/linearalgebra.c \
    -o couple_pendulum -lraylib -lm
```

Or with a simple Makefile target:
```make
CC = gcc
CFLAGS = -O2 -Wall
LIBS = -lraylib -lm

all:
	$(CC) $(CFLAGS) main.c pendulum/pendulum.c integrator/verlet.c \
	    linearalgebra/linearalgebra.c -o couple_pendulum $(LIBS)
```

---

## Running

```bash
./couple_pendulum
```

The program prompts for all parameters interactively:

```
Enter mass of pendulums m1 and m2:          # e.g. 1.0 1.0
Enter length of pendulum l1 and l2:         # e.g. 1.0 1.0
Enter value of spring constant:             # e.g. 2.0
Enter time step:                            # e.g. 0.01
Enter Total time:                           # e.g. 30.0
Enter initial theta of pendulum theta1 and theta2:  # radians, e.g. 0.3 0.0
Enter initial omega of pendulum omega1 and omega2:  # rad/s,   e.g. 0.0 0.0
Enter spring contact point sop1 and sop2:   # 0–1 fraction, e.g. 0.5 0.5
```

Close the window or press **Ctrl+C** to exit.

---

## Parameter Guide

| Parameter | Symbol | Typical Range | Notes |
|-----------|--------|---------------|-------|
| Mass | $m$ | 0.1 – 5.0 kg | Equal masses give symmetric energy transfer |
| Length | $l$ | 0.5 – 2.0 m | Affects natural frequency |
| Spring constant | $k$ | 0.5 – 10.0 N/m | Higher → faster energy exchange |
| Time step | $dt$ | 0.001 – 0.02 s | Smaller is more accurate |
| Spring contact point | $s$ | 0.0 – 1.0 | 1.0 = at bob; 0.5 = midpoint |

---

## Interesting Configurations to Try

**Normal modes** — both pendulums move in phase (or exactly out of phase):
```
theta1=0.3  theta2=0.3   omega1=0  omega2=0    # in-phase
theta1=0.3  theta2=-0.3  omega1=0  omega2=0    # out-of-phase
```

**Beating / energy transfer** — start only one pendulum displaced:
```
theta1=0.3  theta2=0.0   omega1=0  omega2=0
```
Energy will periodically transfer entirely from pendulum 1 to pendulum 2 and back.

---

## License

This project is released for educational use. Feel free to modify and distribute.