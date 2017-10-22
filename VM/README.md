<script type="text/javascript" src="http://cdn.mathjax.org/mathjax/latest/MathJax.js?config=TeX-AMS-MML_HTMLorMML"> </script>

# Vicsek model with quenched disorder

## Model definition
To mimic the collective behaviors of active agents in complex environments, we extend the Vicsek model by dividing the (homogeneous) simulation domain into many square cells with linear size *l*, each of which is assigned a disorder \\(\eta \xi\\), where \\(\epsilon\\) is the strength of the quenched disorder and \\(\xi\\) is a random number uniformly distributed from -PI to PI. The updating rules are:

$$ \begin{align}
\theta^{t+1}_j &= \underbrace{\arg \left[\sum _{k\sim j}e^{i\theta^t_k}\right]} _{\text{alignement}}+\overbrace{\eta\xi^t_j}^{\text{noise}}+\underbrace{\epsilon\xi_m} _{\text{disorder}}, \\\\
\mathbf{x}^{t+1}_j &= \mathbf{x}^t_j+v_0 e^{i\theta^{t+1}_k},
\end{align}$$
where *m* is the index of the cell in which the agent *j* is located at time *t*.

## Usage

## Compiling
Just use `make` command.
```
$ make
```
Then we should get `a.out`. Try `./a.out`, the output is:
```
usage: ./a.out --eta=double --Lx=double --nstep=int [options] ... 
options:
      --eta          noise strength (double)
      --eps          disorder strength (double [=0])
  -L, --Lx           system length in x direction (double)
      --Ly           system length in y direction (double [=0])
  -n, --nstep        total steps to run (int)
  -s, --seed         seed of random number (unsigned long long [=1])
      --rho0         density (double [=0])
  -f, --file         input file (string [=])
      --idx_frame    which frame to read (int [=-1])
      --snap_mode    mode to output snapshots (string [=one])
      --snap_dt      time interval to output snap (int [=2000])
      --log_dt       time interval to log (int [=100000])
      --phi_dt       time interval to calculate phi (int [=100])
  -i, --ini_mode     initializing mode (string [=rand])
      --t_equil      time to reach equilibrium state (int [=100000])
      --cg_dt        time interval to coarse grain (int [=0])
      --cg_ncol      number of cols for coarse grain (int [=0])
      --cg_nrow      number of rows for coarse grain (int [=0])
      --cg_lx        Box size in x for coarse grain (double [=0])
      --cg_ly        Box size in y for coarse grain (double [=0])
      --cg_exp       the exponent for generating frames in log scales (double [=0])
      --cg_format    file format for coarse grain (string [=Hff])
      --cg_win       generate frames in block
      --Sk           calculate structure factor
      --lBox         Box size to calculate correlation (double [=0])
  -?, --help         print this message
```

