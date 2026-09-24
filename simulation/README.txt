This folder contains MATLAB code used to generate panels B, C, and D of
Figure 1 in the associated manuscript.

Files:
- BMC_FIG_1.m   : single script that generates the epidemic simulations
                  shown in Figure 1 (panels B, C, D). Includes a local
                  function, SEIR_ERLANG_ODE_v2, which defines the ODE
                  system for the SEIR model with Erlang-distributed
                  E and I compartments (used with MATLAB's ODE solvers,
                  e.g. ode45).

Requirements:
- MATLAB [version, e.g. R2022a or later]

To run:
Open BMC_FIG_1.m in MATLAB and run. Output is a set of figures
corresponding to Figure 1, panels B-D: a histogram with fitted
Exponential/Erlang distributions (B), and log daily incidence and
cumulative incidence trajectories comparing the SEIR and SE4IR models
under a mid-epidemic transmission reduction (C, D).