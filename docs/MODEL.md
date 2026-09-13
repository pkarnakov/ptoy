# The physical model

What the simulation computes, and why the time integrator and the contact
damping are the way they are. Constants named here live at the top of
`src/particles.cpp`.

## Particles and forces

Every particle has a position, a velocity and a mass `kMass`. The forces acting
on one are summed in `calc_forces()`:

- **Pair force** with every particle closer than $2 r_0$, a repulsive spring
  plus a contact dashpot, from `F12()`.
- **Gravity** $m\mathbf{g}$, when `gravity_enable_` is set.
- **Linear drag** $-\lambda m \mathbf{v}$, with $\lambda$ = `kDissipation`.
- **Mouse force**, attractive or repulsive, when the pointer is held down.
- **Walls**, as four `line` environment objects following the domain.
- **Bonds**, **frozen particles** and **portals**, each in their own pass after
  the block loop.

The repulsive part is a Lennard-Jones-like kernel, cut off at the minimum:

```math
F_\text{spring}(r) = \frac{\sigma}{r}
  \left[ \left(\frac{R}{r}\right)^{12} - \left(\frac{R}{r}\right)^{6} \right],
\qquad R = 2 r_0
```

$F_\text{spring}$ is zero at $r = R$ and `std::max(0., ...)` clamps everything
beyond it, so the force is **exactly** zero past the cutoff, not merely small.
Distance guards that skip pairs past $R$ are therefore bit-exact, which is what
lets `ApplyPortalsForces()` cull pairs without changing results. The dashpot
below is tapered so that it preserves this property.

Particles rest touching at $r = R$. The stiffness there is

```math
k = -\left.\frac{dF}{dr}\right|_{r=R} = \frac{6\sigma}{R^2} = 3750
```

taking the derivative from the compressed side, since the clamp makes the force
one-sided at the cutoff. That stiffness sets the fastest timescale in the
system, and everything below is calibrated against it.

| constant | symbol | value | meaning |
|---|---|---|---|
| `kRadius` | $r_0$ | 0.02 | particle radius, cutoff is $R = 2 r_0$ |
| `kSigma` | $\sigma$ | 1 | pair force strength |
| `kDashpot` | $c$ | 100 | contact damping, see below |
| `kMass` | $m$ | 0.04 | $100\,r_0^2$ |
| `kTimeStep` | $h$ | 5e-4 | the time step |
| `kDissipation` | $\lambda$ | 1e-3 | linear drag rate |
| `kVelocityLimit` | | 10 | hard clamp on speed, applied every half step |
| `kGravity` | $g$ | 10 | default gravity magnitude |

The domain is `[-1,-1]` to `[1,1]` at the reference 800x800 window and follows
the window size from there.

## The contact dashpot

Contacts between hard particles need damping or a pile never settles: grains
trade elastic energy back and forth and the heap shivers forever. The damping
has to be **selective**, though. A drag term strong enough to settle a pile also
visibly slows the bulk flow, which is the motion the player came to watch.

`F12()` therefore adds, on top of the spring,

```math
\mathbf{F}_\text{damp}
  = -c \, \phi(r) \, \big[(\mathbf{v}_1 - \mathbf{v}_2)\cdot\hat{\mathbf{n}}\big]
    \, \hat{\mathbf{n}},
\qquad
\hat{\mathbf{n}} = \frac{\mathbf{p}_1 - \mathbf{p}_2}{r}
```

with the overlap taper

```math
\phi(r) = \max\left(0,\; 1 - \frac{r^2}{R^2}\right)
```

Three properties matter:

- **It acts along the line of centres only.** Shear is untouched, so the
  material still flows and pours like a fluid; only the normal approach and
  rebound of a contact is damped.
- **It is zero unless two particles actually overlap**, and it is zero for any
  motion that carries both particles together, since it sees only
  $\mathbf{v}_1 - \mathbf{v}_2$. Bulk translation, rotation and shear all pass
  through unaffected. This is what a linear drag $-\lambda m\mathbf{v}$ cannot
  do.
- **The taper $\phi$ makes it vanish with the contact.** A constant
  coefficient, as in textbook DEM, jumps to zero the instant a contact breaks;
  tapering keeps the total pair force continuous and keeps it exactly zero past
  $r = R$, preserving the bit-exact cutoff the portal code relies on.

Written as a damping ratio for a pair, with reduced mass $\mu = m/2$:

```math
\zeta = \frac{c \, \phi(r)}{2\sqrt{k\mu}},
\qquad \sqrt{k\mu} = 8.66
```

| overlap | $\phi(r)$ | $\zeta$ | $\gamma h$ |
|---|---|---|---|
| $r = 0.99 R$ | 0.020 | 0.12 | 0.018 |
| $r = 0.98 R$ | 0.040 | 0.23 | 0.035 |
| $r = 0.95 R$ | 0.098 | 0.56 | 0.086 |
| $r = 0.90 R$ | 0.190 | 1.10 | 0.168 |

So a resting contact is lightly damped and springy, while a hard impact is
driven to critical damping ($\zeta = 1$) and stops dead. That is the behaviour
you want from a sandbox pile, and it falls out of the taper rather than being
tuned in.

### It is explicit, and that is fine

The dashpot force is evaluated from the same velocities the force loop already
has, so it is fully explicit and costs no solve. The usual worry applies —
explicit damping has its own stability limit — but the numbers are not close:

- Measured stability limit of the scheme: $\gamma h$ up to 0.99.
- What the contacts actually need: $\gamma h$ between 0.018 and 0.17.

That is a margin of 6x even at heavy overlap. The reason is that the *spring*
is the stiff part of a contact, not the damping: $\omega h$ = 0.15, where
$\omega = \sqrt{k/\mu}$ is the contact frequency, while the damping rate is an
order of magnitude smaller. The dashpot also tightens the step limit on the
spring slightly, from $\omega h \le 2$ to 1.81 at $\zeta$ = 0.1 and 1.49 at
$\zeta$ = 0.3 — we run at 0.15, so this never binds.

Raising $c$ far past its current value does eventually break this. Measured
with a settling pile: 100 is stable and settles well, 200 is marginal, 400 and
above explode — consistent with $\gamma h$ approaching 1.

## The time integrator

Writing $a = F/m$, `Particles::step()` performs:

```math
\begin{aligned}
v^{*}   &= v^{n} + \tfrac{h}{2}\, a(x^{n}, v^{n}) \\
x^{*}   &= x^{n} + \tfrac{h}{2}\, v^{*} \\
v^{n+1} &= v^{n} + h \, a(x^{*}, v^{*}) \\
x^{n+1} &= x^{n} + \tfrac{h}{2}\,\left(v^{n} + v^{n+1}\right)
\end{aligned}
```

A midpoint force evaluation with a trapezoidal position update. The first two
lines are the predictor; note that $x^{*}$ advances with $v^{*}$ rather than
$v^{n}$. The corrector restarts from the state saved in `position_tmp` and
`velocity_tmp` rather than continuing from the predictor.

Both force evaluations take a velocity, because of the dashpot, and the
midpoint rule wants every argument at the midpoint: the corrector uses
$a(x^{*}, v^{*})$, not $a(x^{*}, v^{n})$.

There is a trap in how the code expresses that. `calc_forces()` has no velocity
parameter — it reads `data.velocity`, which holds $v^{*}$ only because the
predictor overwrote it in place beforehand. The requirement is invisible at the
call site, so reordering those loops, or computing the corrector force from a
saved copy of the velocity, would quietly make the scheme first order.

### Accuracy and stability

The scheme is **second order in both position and velocity**, including with
the velocity-dependent dashpot force. Measured convergence at fixed $T$ on a
damped oscillator: the error drops by 4x when $h$ halves, giving 2.04 in $x$
and 2.00 in $v$. Degrading either the position average or the dashpot argument
as above gives 1.00, as expected.

The linear stability limit is $\omega h \le 2$ for the undamped spring, reduced
mildly by the dashpot as noted above. `kVelocityLimit` is a further nonlinear
safety net, clamping speed after each of the two velocity updates.

The scheme is not symplectic. On the harmonic test problem $a(x) = -\omega^2 x$,
with $z = \omega h$, one step is a linear map whose amplification matrix in the
scaled coordinates $(x, v/\omega)$ satisfies

```math
\det M = 1 - \frac{z^4}{8}
```

so it still bleeds a little phase-space volume, 0.007 % per step at contact
frequency. That is a rounding error next to the dashpot and is not relied on.

The leak is there because the scheme is not time-symmetric either. The
trapezoidal position update is, but the velocity update sits on an *explicit*
midpoint predictor, which is not. A symmetric method would need
$\det M(h) \det M(-h) = 1$, and ours is even in $h$, so that would force
$z^4 = 0$. Symmetry implies even order but not the reverse: this scheme is
second order without being symmetric, as explicit midpoint RK2 also is. A
genuinely symmetric scheme, such as implicit midpoint, has $\det M = 1$ and
loses nothing.

### History: the scheme used to be first order on purpose

The position update was originally $x^{n+1} = x^{n} + h\,v^{n+1}$ — using the
new velocity rather than the average. That is a one-term Taylor error, and it
gave

```math
\det M = 1 - \frac{z^2}{2} \qquad \text{exactly}
```

a phase-space contraction of 1.2–2.3 % per step at contact frequency against
0.0005 % for bulk motion at a 1 s period. The scheme was accidentally a
stiffness-selective damper, about 12000x stronger at contact scale than the
linear drag, and that numerical damping was the only thing settling piles. It
was documented as load-bearing and not to be "corrected".

It has since been corrected, because the damping is now explicit. Fixing the
position update alone was tried first and is **not** sufficient: the residual
$z^4/8$ damping is 170x too weak, and a pile that used to settle to a kinetic
energy of about 1 instead plateaued near 116 and shivered indefinitely. Only
with the dashpot added does the second-order scheme settle as well as the old
one:

| | old, 1st order | 2nd order alone | 2nd order + dashpot |
|---|---|---|---|
| settled grid, $E_k$ | 1.41 | 13.7 | **1.02** |
| 1 s after a 10x10 block lands | 13.3 | 182 | **7.3** |
| 2 s after | 1.18 | 116, not decaying | **3.1** |

The dashpot costs about 14 % of a step at 6 threads (69 to 79 µs per step on a
45x45 grid), which is the price of the two extra velocity loads in
`CalcForceAvx`.

The old position update is still reachable: set `kTrapezoidPosition` to false in
`src/particles.cpp` and rebuild. It keeps the dashpot and adds the old numerical
damping on top, so piles settle a few seconds sooner at the cost of first-order
accuracy. It is not a cure for an over-packed domain — pack more particles into
the domain than hexagonal packing holds, for instance by shrinking the window
onto a full pile, and both schemes explode alike.

### Reproducing the numbers

The determinant, the order and the dashpot stability limit are all
straightforward to check on a scalar damped oscillator
$\ddot{x} = -\omega^2 x - 2\gamma\dot{x}$:

```c
/* One step of the ptoy scheme. */
void ptoy_step(double w, double g, double h, double *x, double *v) {
  double a1 = -w*w*(*x) - 2*g*(*v);
  double v_tmp = *v, x_tmp = *x;
  double vs = *v + a1*(h*0.5);          /* v* */
  double xs = *x + vs*(h*0.5);          /* x* */
  double a2 = -w*w*xs - 2*g*vs;         /* dashpot reads v*, not v_tmp */
  *v = v_tmp + a2*h;
  *x = x_tmp + (v_tmp + *v)*(h*0.5);    /* trapezoid, not h * (*v) */
}
```

Stepping the basis vectors $(1, 0)$ and $(0, \omega)$ once each gives the
columns of $M$ in the scaled coordinates $(x, v/\omega)$, and with $\gamma = 0$
its determinant reproduces $1 - z^4/8$ to machine precision. The spectral
radius of the same matrix as a function of $\gamma h$ gives the stability
limit. Integrating to a fixed $T$ and halving $h$ gives the convergence rates.

The pile numbers come from a settling test: build a 25x25 grid with gravity on,
drop a 10x10 block onto it after 1 s with `AddParticleBlock()`, and sum
$\tfrac{1}{2}v^2$ over all particles each frame.
