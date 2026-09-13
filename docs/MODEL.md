# The physical model

What the simulation computes, and why the time integrator is the way it is.
Constants named here live at the top of `src/particles.cpp`.

## Particles and forces

Every particle has a position, a velocity and a mass `kMass`. The forces acting
on one are summed in `calc_forces()`:

- **Pair force** with every particle closer than `2 * kRadius`, from `F12()`.
- **Gravity** `gravity_ * kMass`, when `gravity_enable_` is set.
- **Linear drag** `-v * kDissipation * kMass`.
- **Mouse force**, attractive or repulsive, when the pointer is held down.
- **Walls**, as four `line` environment objects following the domain.
- **Bonds**, **frozen particles** and **portals**, each in their own pass after
  the block loop.

The pair force is a Lennard-Jones-like kernel, repulsive only:

```
F(r) = sigma * ((R/r)^12 - (R/r)^6) / r      with R = 2 * kRadius
```

`F` is zero at `r = R` and `std::max(0., ...)` clamps everything beyond it, so
the force is **exactly** zero past the cutoff, not merely small. Distance
guards that skip pairs past `R` are therefore bit-exact, which is what lets
`ApplyPortalsForces()` cull pairs without changing results.

Particles rest touching at `r = R`. The stiffness there is

```
k = -dF/dr |r=R = 6 * sigma / R^2 = 3750
```

taking the derivative from the compressed side, since the clamp makes the force
one-sided at the cutoff. That stiffness sets the fastest timescale in the system and, below, the amount of
numerical damping the integrator applies.

| constant | value | meaning |
|---|---|---|
| `kRadius` | 0.02 | particle radius, cutoff is `2 * kRadius` |
| `kSigma` | 1 | pair force strength |
| `kMass` | 0.04 | `kRadius^2 * 100` |
| `kTimeStep` | 5e-4 | `dt` |
| `kDissipation` | 1e-3 | linear drag rate |
| `kVelocityLimit` | 10 | hard clamp on speed, applied every half step |
| `kGravity` | 10 | default gravity magnitude |

The domain is `[-1,-1]` to `[1,1]` at the reference 800x800 window and follows
the window size from there.

## The time integrator

`Particles::step()` is **not** forward Euler, and it is **not** leapfrog or
velocity Verlet. Writing `a = F/m` and `h = dt`, one step is:

```
v* = v^n + (h/2) a(x^n)      predictor
x* = x^n + (h/2) v*          (uses v*, not v^n)

v^{n+1} = v^n + h a(x*)      corrector, restarted from velocity_tmp
x^{n+1} = x^n + h v^{n+1}    (uses the new velocity, not v* and not v^n)
```

That is a **midpoint force evaluation with symplectic-Euler-style position
updates**. The corrector restarts from the state saved in `position_tmp` and
`velocity_tmp` rather than continuing from the predictor. Proper midpoint RK2
would advance the position with the midpoint velocity `v*`; using `v^{n+1}`
instead is the one substitution that gives the scheme its character.

### It dissipates energy, by construction

On the harmonic test problem `a(x) = -w^2 x`, with `z = w*h`, the step is a
linear map. In the scaled coordinates `(x, v/w)`, which put both components in
the same units, its amplification matrix is

```
M = [ 1 - z^2 + z^4/4     z (1 - z^2/2) ]
    [ -z (1 - z^2/4)      1 - z^2/2     ]

det M = 1 - z^2 / 2       exactly
```

(The scaling is a diagonal similarity transform, so the determinant is the same
in the raw `(x, v)` coordinates.)

A symplectic integrator has `det M = 1` — velocity Verlet does, which is why
its energy oscillates but does not drift. This scheme contracts phase-space
volume by `1 - (w*h)^2/2` every step. **The damping is in the integrator, not
in the forces.**

### The damping is stiffness-selective, which is why it works here

The loss per step goes as `w^2`, so it is not a uniform viscosity. With the
constants above:

| mode | `w` | `z = w*h` | energy loss per step |
|---|---|---|---|
| particle contact | 306–433 | 0.15–0.22 | **1.2–2.3 %** |
| bulk motion, ~1 s period | 6.3 | 0.003 | 0.0005 % |

(The contact range spans taking the mass as `kMass` or as the reduced mass of a
pair, `kMass/2`.)

Contact vibrations e-fold in about 85 steps, roughly 0.04 game-seconds, while
bulk flow is left almost untouched. That is exactly the right bias for a
contact-dominated system: the stiff modes are the ones that blow up, and they
are the ones being crushed, so the fluid stays stable without looking viscous.
A drag term large enough to do the same job would have damped the visible
motion too.

For comparison, the explicit `kDissipation` term costs `2 * kDissipation * h =
1e-6` of energy per step. At contact scale the integrator damps about
**12000 times harder** than the constant that exists for the purpose.
`kDissipation` matters only for slow bulk motion, where it is within an order
of magnitude of the integrator (0.2 % versus 1 % per game-second).

### Accuracy and stability

Taking one step from exact data, the local truncation errors are uneven:

- **Position: `O(h^2)`**, i.e. first-order consistent. `x^{n+1}` picks up
  `h^2 a` where the true expansion has `h^2 a / 2` — a systematic overshoot
  along the force. Under constant gravity this is a bounded velocity offset of
  `h*a/2`, about 0.0025 at the default settings, i.e. negligible.
- **Velocity: `O(h^3)`**, i.e. second-order consistent. Because `x*` is a
  genuine midpoint estimate, `v^{n+1} = v + h a + (h^2/2) v·grad a + O(h^3)`
  matches the Taylor series one term further than the position does.

**The scheme is globally first order in both variables regardless.** The extra
order in the velocity update does not survive, because the position it is fed
carries `O(h)` error into `a(x*)`. Measured convergence at fixed `T`: the error
in both `x` and `v` halves when `h` halves. The asymmetry is only visible in
the local error — it is a property of the update, not of the solution.

The linear stability limit is `w*h <= 2`, the same as velocity Verlet. The
difference is that Verlet is only marginally stable up to that bound while this
scheme is bleeding energy the whole way, so jammed or overlapping particles
recover instead of heating up. `kVelocityLimit` is a further nonlinear safety
net, clamping speed after each of the two velocity updates.

### If you ever want the damping as a knob

Replacing the predictor coefficient `dt * 0.5` with `dt * theta` in both
predictor lines generalizes the determinant to

```
det M = 1 - theta * z^2
```

so `theta = 0` recovers a symplectic (undamped) Euler-like limit, and the
current `theta = 1/2` sits mid-range. Verified numerically for
`theta` in {0, 1/4, 1/2, 3/4, 1} and `z` in {0.1, 0.5, 1}.

### Do not "fix" this back to Verlet

The scheme looks like a half-finished Verlet and invites being corrected into
one. Doing so removes the numerical damping that is currently holding the
simulation together, and the stability at contact would have to be bought back
with an explicit viscosity that also slows down the motion the player sees.
If the integrator is ever changed, check the stiff-contact behaviour — drop a
block onto a packed grid with `d` and watch whether the pile settles or
shivers.

### Reproducing the numbers

The determinant and the energy decay are straightforward to check by applying
the scheme to a scalar oscillator:

```c
/* One step of the ptoy scheme for a(x) = -w^2 x. */
void ptoy_step(double w, double h, double *x, double *v) {
  double a1 = -w*w*(*x);
  double v_tmp = *v, x_tmp = *x;
  double vs = *v + a1*(h*0.5);   /* v* */
  double xs = *x + vs*(h*0.5);   /* x* */
  double a2 = -w*w*xs;
  *v = v_tmp + a2*h;
  *x = x_tmp + (*v)*h;
}
```

Stepping the basis vectors `(1,0)` and `(0,w)` once each gives the columns of
`M` above in scaled coordinates, and its determinant reproduces `1 - z^2/2` to
machine precision; the same procedure on velocity Verlet gives 1. Integrating
to a fixed `T` and halving `h` gives the convergence rates, and one step from
exact data gives the local truncation errors.
