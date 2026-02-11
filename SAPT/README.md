# Calculation of the LJ scaling parameter

Simply adding explicit induced dipoles to an existing fixed-charge force field will lead to double-counting the induction terms, since these were implicit in the empirical Lennard-Jones potentials.

SAPT calculations can be used to resolve the induction and dispersion terms, and thus calculate a scaling factor to apply to the LJ potential, so that only dispersion contributions remain. The induction energy will be represented explicitly in the polarizable force field.

In most cases, a much less expensive estimation of the scaling factor can be performed [[10.1021/acs.jctc.9b00689]](https://doi.org/10.1021/acs.jctc.9b00689) using only simple quantities of molecular fragments (charge, dipole, polarizability). But in specific cases a full SAPT calculation can be justified.

## Contents

The molecular or ionic fragments in the example provided are $\mathrm{C_2C_1im^+}$ and $\mathrm{C_4H_{10}}$.

SAPT calcualtions are performed using the [Psi4](https://psicode.org) quantum chemistry code, which can be easily installed with conda. Choice of theoretical level and basis set followed the work of Sherrill [10.1063/1.4867135](https://doi.org/10.1063/1.4867135)

* `c2c1im_c4_opt.in`: geometry optimization of the pair of fragments
* `c2c1im_c4_sapt0.in`: faster calculations at the SAPT0 level to obtain an approximate potential energy surface.
* `c2c1im_c4_sapt2.in`: full second-order SAPT2+ energies, to perform near the minimum.
* `sapt_out.sh`: script to extract the potential energy terms from Psi4 output.

In the input files above, the total energy is printed with a tilde '~', enabling the easy use of `grep` to retrive the value. 
