# MCSMD
Multi-centred semiclassical Moyal dynamics

This is a developing method:
* Version 0: basic least square fit with prespecified centre(s)
* Version 1: add normalization and purity conservation constraint
* Version 2: automatic centre generation and merge

Supported job types:
1. SMD
* Propagate SMD quantities based on user specified initial condition, and save the SMD quantities and Wigner coefficients
2. Wigner
* Compute Wigner distribution to plot based on saved Wigner coefficients

Analyze.py provides visualization for each type of job. Just run it directly after MCSMD.exe finishes

Reference:
> 1. Y. Shen, L. Wang 2018 J. Chem. Phys.