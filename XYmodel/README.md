<script type="text/javascript" src="http://cdn.mathjax.org/mathjax/latest/MathJax.js?config=TeX-AMS-MML_HTMLorMML"> </script>

# XY model

In the XY model, the phase \\(\phi_i\\) at each lattice site *i* evolves according to
$$
\gamma\frac{d\phi_i}{dt}=\kappa\sum_j \sin(\phi_i-\phi_j)-\eta_i,
$$
where \\(\eta_i\\) is the Langevin noise at site *i*, which has zero mean and the moments
$$
\langle\eta_i(t)\eta_j(t')\rangle=2k_B \gamma \delta_{i,j} \delta (t-t').
$$
The numerical integration procedure is Euler update
$$
\phi_i(t+\Delta t) = \phi_i(t)-\Delta t \left [\eta_i(t)+\frac{\kappa}{\gamma}\sum_j \sin\left (\phi_i(t)-\phi_j(t)\right )\right ].
$$

The simulations are performed on a square lattice and the sum over *j* was carried out over the eight nearest neighbors of *i*. The constans \\(\kappa\\) and \\(\gamma\\) can be scaled and are thus taken to be 1. The time step is taken to be 0.05. The Langevin noise \\(\eta_i(t)=2\pi c_L r_i\\) consisted of random numbers \\(r_i\\) having a uniform distribution over the interval from -0.5 to 0.5. The order-to-disorder transition is determined to occur near \\(c_L\sim3\\).

Ref:

B. Yurke. et al., *Coarsening dynamics of the XY model.* 1993 PRE