# The physical model

What the simulation computes, and why the time integrator and the contact
damping are the way they are. Constants named here live at the top of
`src/particles.cpp`.

## Particles and forces

Every particle has a position, a velocity and a mass `kMass`. The forces acting
on one are summed in `calc_forces()`:

- **Pair force** with every particle closer than `2 * kRadius`, a repulsive
  spring plus a contact dashpot, from `F12()`.
- **Gravity** `gravity_ * kMass`, when `gravity_enable_` is set.
- **Linear drag** `-v * kDissipation * kMass`.
- **Mouse force**, attractive or repulsive, when the pointer is held down.
- **Walls**, as four `line` environment objects following the domain.
- **Bonds**, **frozen particles** and **portals**, each in their own pass after
  the block loop.

The repulsive part is a Lennard-Jones-like kernel, cut off at the minimum:

```
F_spring(r) = sigma * ((R/r)^12 - (R/r)^6) / r      with R = 2 * kRadius
```

`F_spring` is zero at `r = R` and `std::max(0., ...)` clamps everything beyond
it, so the force is **exactly** zero past the cutoff, not merely small.
Distance guards that skip pairs past `R` are therefore bit-exact, which is what
lets `ApplyPortalsForces()` cull pairs without changing results. The dashpot
below is tapered so that it preserves this property.

Particles rest touching at `r = R`. The stiffness there is

```
k = -dF/dr |r=R = 6 * sigma / R^2 = 3750
```

taking the derivative from the compressed side, since the clamp makes the force
one-sided at the cutoff. That stiffness sets the fastest timescale in the
system, and everything below is calibrated against it.

| constant | value | meaning |
|---|---|---|
| `kRadius` | 0.02 | particle radius, cutoff is `2 * kRadius` |
| `kSigma` | 1 | pair force strength |
| `kDashpot` | 100 | contact damping, see below |
| `kMass` | 0.04 | `kRadius^2 * 100` |
| `kTimeStep` | 5e-4 | `dt` |
| `kDissipation` | 1e-3 | linear drag rate |
| `kVelocityLimit` | 10 | hard clamp on speed, applied every half step |
| `kGravity` | 10 | default gravity magnitude |

The domain is `[-1,-1]` to `[1,1]` at the reference 800x800 window and follows
the window size from there.

## The contact dashpot

Contacts between hard particles need damping or a pile never settles: grains
trade elastic energy back and forth and the heap shivers forever. The damping
has to be **selective**, though. A drag term strong enough to settle a pile also
visibly slows the bulk flow, which is the motion the player came to watch.

`F12()` therefore adds, on top of the spring,

```
F_damp = -kDashpot * overlap(r) * ((v1 - v2) . n) n        n = (p1 - p2) / r

overlap(r) = max(0, 1 - r^2 / R^2)
```

Three properties matter:

- **It acts along the line of centres only.** Shear is untouched, so the
  material still flows and pours like a fluid; only the normal approach and
  rebound of a contact is damped.
- **It is zero unless two particles actually overlap**, and it is zero for any
  motion that carries both particles together. Bulk translation, rotation and
  shear all pass through it unaffected. This is what `kDissipation` cannot do.
- **The `overlap` taper makes it vanish with the contact.** A constant
  coefficient, as in textbook DEM, jumps to zero the instant a contact breaks;
  tapering keeps the total pair force continuous and keeps it exactly zero past
  `r = R`, preserving the bit-exact cutoff the portal code relies on.

Written as a damping ratio for a pair, with reduced mass `mu = kMass / 2`:

```
zeta = kDashpot * overlap / (2 * sqrt(k * mu))        sqrt(k * mu) = 8.66
```

| overlap | `overlap(r)` | `zeta` | `gamma * h` |
|---|---|---|---|
| `r = 0.99 R` | 0.020 | 0.12 | 0.018 |
| `r = 0.98 R` | 0.040 | 0.23 | 0.035 |
| `r = 0.95 R` | 0.098 | 0.56 | 0.086 |
| `r = 0.90 R` | 0.190 | 1.10 | 0.168 |

So a resting contact is lightly damped and springy, while a hard impact is
driven to critical damping and stops dead. That is the behaviour you want from
a sandbox pile, and it falls out of the taper rather than being tuned in.

### It is explicit, and that is fine

The dashpot force is evaluated from the same velocities the force loop already
has, so it is fully explicit and costs no solve. The usual worry applies —
explicit damping has its own stability limit — but the numbers are not close:

- Measured stability limit of the scheme in `gamma * h`: **0.99**.
- What the contacts actually need: `gamma * h` = **0.018 to 0.17**.

That is a margin of 6x even at heavy overlap. The reason is that the *spring* is
the stiff part of a contact, not the damping: `w * h = 0.15` where `w` is the
contact frequency, while the damping rate is an order of magnitude smaller. The
dashpot also tightens the step limit on the spring slightly, from `w*h <= 2.0`
to `w*h <= 1.81` at `zeta = 0.1` and `1.49` at `zeta = 0.3` — we run at 0.15,
so this never binds.

Raising `kDashpot` far past its current value does eventually break this.
Measured with a settling pile: 100 is stable and settles well, 200 is marginal,
400 and above explode — consistent with `gamma * h` approaching 1.

## The time integrator

Writing `a = F/m` and `h = dt`, `Particles::step()` performs:

```
v* = v^n + (h/2) a(x^n, v^n)      predictor
x* = x^n + (h/2) v*               (uses v*, not v^n)

v^{n+1} = v^n + h a(x*, v*)       corrector, restarted from velocity_tmp
x^{n+1} = x^n + (h/2) (v^n + v^{n+1})
```

A midpoint force evaluation with a trapezoidal position update. The corrector
restarts from the state saved in `position_tmp` and `velocity_tmp` rather than
continuing from the predictor.

Two details are load-bearing and easy to undo by accident:

- **The dashpot in the corrector is evaluated at `v*`, not `v^n`.** The code
  gets this for free because the predictor writes `v*` into `data.velocity`
  in place before `calc_forces()` runs again. Feeding `v^n` there instead drops
  the scheme to first order (measured: 1.00 instead of 2.00).
- **The position update averages the two velocities.** Using `v^{n+1}` alone
  overshoots by `h^2 a / 2` and is also first order.

### Accuracy and stability

The scheme is **second order in both position and velocity**, including with
the velocity-dependent dashpot force. Measured convergence at fixed `T` on a
damped oscillator: the error drops by 4x when `h` halves.

| position update | dashpot velocity | global order `x`, `v` |
|---|---|---|
| `x + h v^{n+1}` | `v^n` | 1.00, 0.99 |
| `x + h v^{n+1}` | `v*` | 1.00, 1.00 |
| `x + (h/2)(v^n + v^{n+1})` | `v^n` | 1.00, 1.00 |
| `x + (h/2)(v^n + v^{n+1})` | **`v*`** | **2.04, 2.00** |

The linear stability limit is `w*h <= 2` for the undamped spring, reduced
mildly by the dashpot as noted above. `kVelocityLimit` is a further nonlinear
safety net, clamping speed after each of the two velocity updates.

The scheme is not symplectic: on a harmonic oscillator its amplification matrix
has

```
det M = 1 - z^4 / 8        z = w * h
```

so it still bleeds a little phase-space volume, 0.007 % per step at contact
frequency. That is a rounding error next to the dashpot and is not relied on.

### History: the scheme used to be first order on purpose

The position update was originally `x^{n+1} = x^n + h v^{n+1}` — using the new
velocity rather than the average. That is a one-term Taylor error, and it gave

```
det M = 1 - z^2 / 2        exactly
```

a phase-space contraction of **1.2–2.3 % per step** at contact frequency
against 0.0005 % for bulk motion at a ~1 s period. The scheme was accidentally
a stiffness-selective damper, about 12000x stronger at contact scale than
`kDissipation`, and that numerical damping was the only thing settling piles.
It was documented as load-bearing and not to be "corrected".

It has since been corrected, because the damping is now explicit. Fixing the
position update alone was tried first and is **not** sufficient: the residual
`z^4 / 8` damping is 170x too weak, and a pile that used to settle to a kinetic
energy of ~1 instead plateaued at ~116 and shivered indefinitely. Only with the
dashpot added does the second-order scheme settle as well as the old one:

| | old, 1st order | 2nd order alone | 2nd order + dashpot |
|---|---|---|---|
| settled grid, KE | 1.41 | 13.7 | **1.02** |
| 1 s after a 10x10 block lands | 13.3 | 182 | **7.3** |
| 2 s after | 1.18 | 116, not decaying | **3.1** |

The dashpot costs about 14 % of a step at 6 threads (69 -> 79 us/step on a
45x45 grid), which is the price of the two extra velocity loads in
`CalcForceAvx`.

### Reproducing the numbers

The determinant, the order and the dashpot stability limit are all
straightforward to check on a scalar damped oscillator `x'' = -w^2 x - 2 g x'`:

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

Stepping the basis vectors `(1,0)` and `(0,w)` once each gives the columns of
the amplification matrix in the scaled coordinates `(x, v/w)`, and with `g = 0`
its determinant reproduces `1 - z^4/8` to machine precision. The spectral
radius of the same matrix as a function of `g*h` gives the stability limit.
Integrating to a fixed `T` and halving `h` gives the convergence rates.

The pile numbers come from a settling test: build a 25x25 grid with gravity on,
drop a 10x10 block onto it after 1 s with `AddParticleBlock()`, and sum
`0.5 * v^2` over all particles each frame.
