# pol-il
Polarisable force field for ionic liquids
'fftool 1 c4c1im.zmat 1 dca.zmat -b 20
'packmol <pack.inp'
'fftool 1 c4c1im.zmat 1 dca.zmat -b 20 -a -l'

'python polarizer.py c4c1im.zmat dca.zmat -f alpha.ff -q -id data.lmp -od data-p.lmp'

> # Commands to include in the LAMMPS input script
> 
> # adapt the pair_style command as needed
> pair_style hybrid/overlay [...] coul/long/cs 12.0 thole 2.600 12.0
> 
> # data file with Drude oscillators added
> read_data data-p.lmp
> 
> # pair interactions with Drude particles written to file
> # Thole damping recommended if more than 1 Drude per molecule
> include pair-drude.lmp
> 
> # atom groups convenient for thermostats (see package documentation), etc.
> group ATOMS type 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15
> group CORES type 1 2 3 4 6 9 10 12 13 14 15
> group DRUDES type 16 17 18 19 20 21 22 23 24 25 26
> 
> # flag for each atom type: [C]ore, [D]rude, [N]on-polarizable
> fix DRUDE all drude C C C C N C N N C C N C C C C D D D D D D D D D D D
> 
> # store velocity information of ghost atoms
> comm_modify vel yes
> 
> # compute the temperatures of ATOMS, DC-DP pair centers of mass and DPs
> compute TATOM ATOMS temp
> compute TDRUDE all temp/drude
> 
> # examples of termostats and barotats
> # NVT (Nose-Hoover)
> fix DTDIR all drude/transform/direct
> fix TSTAT ATOMS nvt temp 300.0 300.0 200
> fix TSTDR DRUDES nvt temp 1.0 1.0 50
> fix DTINV all drude/transform/inverse
> # NPT (Nose-Hoover)
> fix DTDIR all drude/transform/direct
> fix TSTAT ATOMS npt temp 300.0 300.0 200 iso 1.0 1.0 1000
> fix_modify TSTAT temp TATOM press thermo_press
> fix TSTDR DRUDES nvt temp 1.0 1.0 50
> fix DTINV all drude/transform/inverse
> 
> # avoiding the flying ice cube artefact
> fix ICECUBE all momentum 1000 linear 1 1 1
> 
> # output the temperatures of ATOMS, DC-DP pair centers of mass and DPs
> thermo_style custom step [...] c_TATOM c_TDRUDE[1] c_TDRUDE[2]
> 
> # write Drude particles to dump file
> dump_modify ... element ... D D D D D D D D D D D
> 
> # ATTENTION!
> #  * read_data may need 'extra/special/per/atom' keyword, LAMMPS will exit with a message.
> #  * If using fix shake the group-ID must not include Drude particles. Use group ATOMS, for example.
> #  * Give all I<=J pair interactions, no mixing.
> #  * Pair style coul/long/cs from CORESHELL package is used for interactions
> #    of Drude particles. Alternatively pair lj/cut/thole/long could be used,
> #    avoiding hybrid/overlay and allowing mixing. See doc pages.

'cat pair.lmp pair-drude.lmp >  pair-p.lmp'

'python scaleLJ.py -f fragment.ff -a alpha.ff -i fragment.inp -ip pair-p.lmp -op pair-p-sc.lmp'

> Epsilon LJ parameters was scaled by k_pred parameter. Changes are marked with '~'.
> Sigma LJ parameters were not scaled.
> --------------------------------------------------------------------
>  Fragment1	 Fragment2	  k_sapt	 k_pred
>     c2c1im	     c4h10	    0.76 	   0.78
>     c2c1im	       dca	    0.61 	   0.68
>      c4h10	     c4h10	    0.94 	   1.00
>      c4h10	       dca	    0.69 	   0.72
> --------------------------------------------------------------------
