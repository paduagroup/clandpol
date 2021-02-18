#!/usr/bin/env python
# scaleLJ.py - scale epsilon and sigma LJ parameters in pair.lmp LAMMPS file.
# Agilio Padua <agilio.padua@ens-lyon.fr>
# Kateryna Goloviznina <kateryna.goloviznina@ens-lyon.fr>
# version 2021/02/18

import sys
import math
import argparse
import os
import numpy as np

usage = """
==============================================================================
Scale epsilon and sigma LJ parameters in pair.lmp LAMMPS file
------------------------------------------------------------------------------
Format of file containing specification of Drude oscillators 
and polarisability values (alpha.ff):
  # type  dm/u   dq/e  k/(kJ/molA2)  alpha/A3  a_thole
  CT      0.4    -1.0     4184.0      1.016     2.6 
  ...
* dm is the mass to place on the Drude particle (taken from its core),
* dq is the charge to place on the Drude particle (taken from its core),
* k is the harmonic force constant of the bond between core and Drude,
* alpha is the polarizability, hyrdogen aroms are not merged,
* a_thole is a parameter of the Thole damping function.
------------------------------------------------------------------------------
Format of file containing monomers and dimers specification (fragment.ff)
  MONOMERS
  # name       q/e       mu/D
  c2c1im       1.0      1.1558
  ..
  DIMERS
  # m1         m2       r_COM/A    k_sapt
  c2c1im       dca       2.935      0.61
  ...
* q is the charge of the monomer,
* mu is the dipole moment of the monomer,
* m1 and m2 are the monomers forming a dimer,
* r_COM distance between the center of mass of the monomers,
* k_sapt is the scaling factor for the epsilon of LJ potential, 
obtained by SAPT quantum calculation (optional).
------------------------------------------------------------------------------
Format of file containing fragments list with atomic indices (fragment.inp)
  #name  indices
  c2c1im 1:10
  dca    11:13
  ...
* atomic indices or/and a range of indices that correspond to atomic types 
associating with this fragment in data.lmp file. 
------------------------------------------------------------------------------
Script requires structure files of monomers used in fragment.inp file 
(.zmat, .xyz, .mol, .pdb)
==============================================================================
"""

class _Const(object):
    
    sigma_k = 0.985
    
    # C0 and C1 are coefficients to predict k
    @staticmethod
    def C0():
        return 0.254952   
    @staticmethod
    def C1():
        return 0.106906

    @staticmethod
    def isfloat(value):
        try:
            float(value)
            return True
        except ValueError:
            return False
    

# monomer in fragment.ff        
class Monomer (object):
    def __init__(self, name, q, mu):
        self.name = name
        self.q = q
        self.mu = mu

    def __str__(self):
        return '%10s  q = %2d   mu = %6.4f D' % (self.name, self.q, self.mu)

# dimer in fragment.ff
class Dimer (object):
    def __init__(self, m1, m2):
        self.m1 = m1
        self.m2 = m2
        self.r = None
        self.k_sapt = None

    def SetR(self, r):
        self.r = r

    def SetKSAPT(self,  k_sapt):
        self.k_sapt = k_sapt 

    def __str__(self):
        res = '%10s %10s' % (self.m1.name, self.m2.name)
        if self.r is not None:
            res += '    r = %6.4f' % self.r
        if self.k_sapt is not None:
            res += '    k_sapt = %4.2f' % self.k_sapt
        return res

class Forcefield(object):
    #  read fragment.ff
    def __init__(self, filename):
        self.filename = filename
        self.monomers = []
        self.dimers= []
        try:
            with open(self.filename, 'r') as f:
                for line in f:
                    if line.startswith('#') or line.strip() == '':
                        continue
                    if line.lower().startswith('monomer'):
                        section = 'monomers'
                        continue
                    elif line.lower().startswith('dimer'):
                        section = 'dimers'
                        continue

                    tok = line.strip().split()

                    if section == 'monomers':
                        if tok[0].endswith('+') or tok[0].endswith('-'):    
                            tok[0] = tok[0][:-1]
                        name = tok[0].lower()
                        q = float(tok[1])
                        mu  = float (tok[2])

                        if next((x for x in self.monomers if x.name == name), None) is not None:
                            raise Exception('  error: monomer ' + name + ' is specified twice in ' + self.filename)

                        self.monomers.append(Monomer(name, q, mu))
                    
                    elif section == 'dimers':
                        for i in range(0,2):
                            if tok[i].endswith('+') or tok[i].endswith('-'):    
                                tok[i] = tok[i][:-1]
                        m1_name = tok[0].lower()
                        m2_name = tok[1].lower()
                        r = float (tok[2])

                        if next((x for x in self.dimers if (x.m1.name == m1_name and x.m2.name == m2_name) or (x.m1.name == m2_name and x.m2.name  == m1_name)), None):
                            raise Exception('  error: dimer ' + m1_name+' '+m2_name + ' is specified twice in ' + self.filename)

                        d = self.SetDimer(m1_name, m2_name)
                        self.dimers.append(d)
                        d.SetR(r)
                        if len(tok)>3:
                            k_sapt = float (tok[3])
                            d.SetKSAPT(k_sapt)

        except IOError:
            print('  error: force field file ' + self.filename + ' not found')
            sys.exit(1)
        except IndexError:
            print('  error: incorrect force field file line:' + line)
            sys.exit(1)
        except Exception as e:
            print(e)
            sys.exit(1)

    # create new dimer from two monomers
    def SetDimer(self, m1_name, m2_name):
        try:
            m1 = next((x for x in self.monomers if x.name == m1_name), None)
            m2 = next((x for x in self.monomers if x.name == m2_name), None)
            if (m1 is not None and m2 is not None):
                d = Dimer(m1, m2)
                return d
            else:
                raise Exception('  error: monomer '+m1_name+' or monomer '+m2_name+' not descibed in monomers section of '+ self.filename)

        except Exception as e:
            print(e)
            sys.exit(1)

    def __str__(self):
        res = self.filename
        res += '\nMONOMERS'
        for m in self.monomers:
            res+='\n'+str(m)
        res += '\nDIMERS'
        for d in self.dimers:
            res+='\n'+str(d)
        return res    

# fragment in fragment.inp; fragment is based on monomer with index range in currect system added        
class Fragment(Monomer):
    def __init__(self, m, ind_range, res):
        self.name = m.name
        self.q = m.q
        self.mu = m.mu
        self.ind_range = ind_range
        self.pol_model = res[0]
        self.scale_eps = res[1]
        self.scale_sig = False
    
    def __str__(self):
            return '%10s  q = %5.2f  mu = %6.4f D    Pol = %5s ScaleEps = %5s ScaleSig = %5s Atoms=%s ' % (self.name, self.q, self.mu,self.pol_model, self.scale_eps, self.scale_sig, self.ind_range)
    
    # checks if fragment is polarisable
    @staticmethod
    def PolExclude(m, p):
        pol_model = False
        scale_eps = True

        if p is not None:
            if m.name in p:
                pol_model = True
                if (m.q == 0 and m.mu == 0):
                    scale_eps = False
        
        return (pol_model,scale_eps)

    # read atoms from .mol for a given fragment
    def GetAtomsFromMol(self, pol):
        zfilename = self.name+'.mol'
        with open(zfilename, 'r') as f:
            atoms = []
            tok = f.readline().strip().split()
            self.name = tok[0]
            line = f.readline() 
            line = f.readline()
            tok = f.readline().strip().split()   # counts line
            natom = int(tok[0])
            
            for i in range(natom):
                tok = f.readline().strip().split()
                at_name = tok[3]
                print(at_name)
                atom = next((x for x in pol.atomtypes if x.name == at_name), None)
                if atom is None:
                    raise Exception('  error: atom type '+ at_name + ' not found in '+ pol.filename)
                atoms.append(atom)
        self.atoms = atoms

    # read atoms from .pdb for a given fragment
    def GetAtomsFromPdb(self, pol):
        zfilename = self.name+'.pdb'
        with open(zfilename, 'r') as f:
            atoms = []
            line = f.readline()
            while 'COMPND' not in line:
                line = f.readline()
            tok = line.strip().split()
            self.name = tok[1]
            line = f.readline()
            while not (line.startswith('ATOM') or line.startswith('HETATM')):
                line = f.readline()
            while 'ATOM' in line or 'HETATM' in line:
                tok = line.strip().split()
                at_name = tok[2]
                atom = next((x for x in pol.atomtypes if x.name == at_name), None)
                if atom is None:
                    raise Exception('  error: atom type '+ at_name + ' not found in '+ pol.filename)
                atoms.append(atom)
                line = f.readline()
        self.atoms = atoms

    # read atoms from .zmat for a given fragment
    def GetAtomsFromZmat(self, pol):
        try:
            zfilename = self.name+'.zmat'
            with open(zfilename, 'r') as f:
                atoms = []
                line = f.readline()
                while line.strip().startswith('#'):
                    line = f.readline()
                fr_name = line.strip()
                line = f.readline()
                while line.strip().startswith('#') or line.strip() == '':
                    line = f.readline()
        
                tok = line.strip().split()
                if len(tok) > 1:   # there can be line numbers
                    shift = 1
                else:
                    shift = 0
                while line:
                    tok = line.strip().split()
                    if len(tok) == 0:
                        break
                    at_name = tok[shift]
                    atom = next((x for x in pol.atomtypes if x.name == at_name), None)
                    if atom is None:
                        raise Exception('  error: atom type '+ at_name + ' not found in '+ pol.filename)
                    atoms.append(atom)
                    line = f.readline()
            self.atoms = atoms
                            
        except IOError:
            print('  error: fragment zmat file ' + filename + ' not found')
            sys.exit(1)
        except Exception as e:
            print(e)
            sys.exit(1)

    # read atoms from .xyz for a given fragment
    def GetAtomsFromXYZ(self, pol):
        try:
            xyzfilename = self.name+'.xyz'
            with open(xyzfilename, 'r') as f:
                atoms = []
                line = f.readline()
                while line.strip().startswith('#'):
                    line = f.readline()
                n_atoms = int(line.strip())
                line = f.readline()
                fr_name = line.strip()
                line = f.readline()
                while line.strip().startswith('#') or line.strip() == '':
                    line = f.readline()
        
                tok = line.strip().split()
                while line:
                    tok = line.strip().split()
                    if len(tok) == 0:
                        break
                    at_name = tok[0]
                    atom = next((x for x in pol.atomtypes if x.name == at_name), None)
                    if atom is None:
                        raise Exception('  error: atom type '+ at_name + ' not found in '+ pol.filename)
                    atoms.append(atom)
                    line = f.readline()
            self.atoms = atoms
                            
        except IOError:
            print('  error: fragment zmat file ' + filename + ' not found')
            sys.exit(1)
        except Exception as e:
            print(e)
            sys.exit(1)

# fragment pair that is based on dimer but uses fragments (not monomers) as units            
class FragmentPair(Dimer):
    def __init__(self, fr1, fr2, r, k_sapt):
        self.fr1 = fr1
        self.fr2  = fr2
        self.r = r
        self.k_sapt = k_sapt

    def __str__(self):
        res = '%10s %10s' % (self.fr1.name, self.fr2.name)
        if self.r is not None:
            res += '    r = %6.4f' % self.r
        if self.k_sapt is not None:
            res += '    k_sapt = %4.2f' % self.k_sapt
        return res

    # predicts k factor for a given fragment pair based on total charge, alpha, dipole moment of fragment and interfragment distance 
    def PredictK(self):
        alpha_fr1 = [sum(x.alpha for x in self.fr1.atoms)][0]
        alpha_fr2 = [sum(x.alpha for x in self.fr2.atoms)][0]
        k_pred = 1.0

        if not self.fr1.pol_model:
            k_pred += _Const.C0()*self.r*self.r*(self.fr2.q*self.fr2.q)/alpha_fr2
            k_pred += _Const.C1()*(self.fr2.mu*self.fr2.mu)/alpha_fr2

        if not self.fr2.pol_model:
            k_pred += _Const.C0()*self.r*self.r*(self.fr1.q*self.fr1.q)/alpha_fr1
            k_pred += _Const.C1()*(self.fr1.mu*self.fr1.mu)/alpha_fr1
        
        k_pred = 1/k_pred
        self.k_pred = k_pred

# atom type with polarisability value       
class AtomType(object):
    def __init__(self, name, alpha):
        self.name = name
        self.alpha = alpha
    
    def __str__(self):
        return  '%s    %6.3f' %  (self.name, self.alpha)

# polarisability values for all atom types; from alpha.ff
class Polarisation(object):
    def __init__(self, filename):
        self.filename = filename
        self.atomtypes = []
        try:
            with open(self.filename, "r") as f:
                for line in f:
                    line = line.strip()
                    if line.startswith('#') or len(line) == 0:
                        continue
                    tok = line.split()
                    a_name = tok[0]
                    a_alpha = float(tok[4])
                    if next((x for x in self.atomtypes if x.name == a_name), None) is not None:
                        raise Exception('  error: atom type ' + a_name + ' is specified twice in ' + self.filename)
                    self.atomtypes.append(AtomType(a_name,a_alpha))

        except IOError:
            print('  error: polarisation file ' + self.filename + ' not found')
            sys.exit(1)
        except Exception as e:
            print(e)
            sys.exit(1)
    
    def __str__(self):
        res = self.filename
        for at in self.atomtypes:
            res += '\n' + str(at)
        return res

# system that consist of fragments from fragment.inp 
class System(object):
    def __init__(self, filename):
        self.fragments = []
        self.fragmentpairs = []
        self.filename = filename
    
    # read fragments from fragment.inp    
    def GetFragments(self, ff, p):
        try:
            with open(self.filename, 'r') as f:
                for line in f:
                    if line.startswith('#') or line.strip() == '':
                        continue
                    else:
                        tok = line.strip().split()
                        m_name = tok[0].lower()
                        if next((x for x in self.fragments if x.name == m_name), None) is not None:
                            raise Exception('  error: fragment ' + m_name + ' is specified twice in ' + self.filename)
                        if len(tok)<2:
                            raise Exception('  error: no index range for ' + m_name + ' fragment in ' + self.filename)
                        ind_range = []
                        for ind in tok[1:]:
                            if ':' in ind:
                                tok_ind = ind.split(':')
                                ind_range.extend(range(int(tok_ind[0]),int(tok_ind[1])+1))
                            else:
                                ind_range.append(int(ind))
                    m = next((x for x in ff.monomers if x.name == m_name), None)
                    if m is None:
                        raise Exception('  error: fragment ' + m_name + ' not found in ' + ff.filename)
                    
                    self.fragments.append(Fragment(m,ind_range,Fragment.PolExclude(m,p)))

        except IOError:
            print('  error: fragment input file ' + filename + ' not found')
            sys.exit(1)
        except Exception as e:
            print(e)
            sys.exit(1)

    # generates fragment pairs based on fragment list; fragment pairs are not formed for ++ and -- combinations
    def GetFragmentPairs(self, ff):
        try:
            i = j = 0
            for i in range(0,len(self.fragments)):
                for j in range(i,len(self.fragments)):
                    if (np.sign(self.fragments[i].q) != np.sign(self.fragments[j].q) or max(abs(self.fragments[i].q),abs(self.fragments[j].q)) == 0) and (self.fragments[i].scale_eps and self.fragments[j].scale_eps) and (not self.fragments[i].pol_model or not self.fragments[j].pol_model):
                        d = next((x for x in ff.dimers if (x.m1.name == self.fragments[i].name and x.m2.name == self.fragments[j].name) or (x.m1.name == self.fragments[j].name and x.m2.name == self.fragments[i].name)), None)
 
                        if d is not None:
                            self.fragmentpairs.append(FragmentPair(self.fragments[i],self.fragments[j],d.r,d.k_sapt))
                        else:
                            raise Exception('  error: dimer ' + self.fragments[i].name + ' ' + self.fragments[j].name + ' not found in ' + ff.filename)
                    j+=1
                i+=1
        except Exception as e:
            print(e)  
            sys.exit(1)

    def ParseScaleSigma(self, scsig):
        
        if len(scsig) > 0:
            if _Const.isfloat(scsig[0]):
                _Const.sigma_k = float(scsig[0])
                scsig.pop(0)

        if len(scsig) > 0:
            for i in scsig:
                fi = next((x for x in self.fragments if (x.name == i)), None)
                if (fi is None):
                    raise Exception('  error: fragment %s specified with -s option not found in fragment.inp ' % i)
                elif fi.pol_model:
                    raise Exception('  error: fragment %s specified with -s option is already polarisable: sigma should not be scaled' % i)
                else:
                    fi.scale_sig = True
        else:
            for f in self.fragments:
                f.scale_sig = True 

   
    def __str__(self):
        res = self.filename
        res+='\nFRAGMENTS'
        for fr in self.fragments:
            res+='\n'+str(fr)
        res+='\nFRAGMENT PAIRS'
        for frp in self.fragmentpairs:
            res+='\n'+str(frp)
        return res

    # call a function to get atom list for fragments based on existing atom type 
    def GetFragAtoms(self, pol):
        try:
            for fr in self.fragments:
                if os.path.exists('./'+fr.name+'.zmat'):
                    fr.GetAtomsFromZmat(pol)
                elif os.path.exists('./'+fr.name+'.xyz'):
                    fr.GetAtomsFromXYZ(pol)
                elif os.path.exists('./'+fr.name+'.mol'):
                    fr.GetAtomsFromMol(pol)
                elif os.path.exists('./'+fr.name+'.pdb'):
                    fr.GetAtomsFromPdb(pol)
                else:
                    raise Exception('  error: structure file (.zmat, .xyz, .pdb or .mol) for fragment '+ fr.name + ' not found')
         
        except Exception as e:
            print(e)
            sys.exit(1)

    def GetKPred(self):
        for frp in self.fragmentpairs:
            frp.PredictK()
 
# class with static functions to scale sigma and epsilon in pair-p.lmp and print output file             
class ScaleLJ(object):
    @staticmethod
    def Scale(pair_in_file, ff, syst, sapt, scsig):
        try:
            res= []
            for line in open(pair_in_file, 'r'):
                tok = line.strip().split()
                if tok[0] != 'pair_coeff':
                    continue
                if not (tok[3].startswith('lj') and tok[1].isdigit() and tok[2].isdigit()):
                    res.append(line.strip())
                    continue
                i = int(tok[1])
                j = int(tok[2])
                pair = tok[3]
                eps = float(tok[4])
                sig = float(tok[5])
                com = ''
                if len(tok) >= 7:
                    for n in range(6, len(tok)):
                        com += ' ' + tok[n]                   
                for frp in syst.fragmentpairs:
                    if not sapt:
                        k = frp.k_pred
                    elif frp.k_sapt is not None:
                        k = frp.k_sapt
                    else:
                        raise Exception('  error: k_sapt for '+ frp.fr1.name +' '+frp.fr2.name + ' dimer not found in '+ ff.filename)                   
                    if ((i in frp.fr1.ind_range and j in frp.fr2.ind_range) or (j in frp.fr1.ind_range and i in frp.fr2.ind_range)) and k < 1:
                        eps *= k
                        com += ' ~'
                        break                    
                
                if (scsig is not None):
                    fi = next((x for x in syst.fragments if (i in x.ind_range)), None)
                    fj = next((x for x in syst.fragments if (j in x.ind_range)), None)

                    if (fi is not None and fj is not None):
                        if fi.scale_sig and fj.scale_sig:
                            sig *= _Const.sigma_k
                            com += ' *'
  
                res.append("pair_coeff {0:4d} {1:4d} {2:18s} {3:10.6f}   {4:10.6f} {5:s}".format(i, j, pair, eps, sig, com))
            return res
        
        except IOError:
            print('  error: pair style file ' + pair_file + ' not found')
            sys.exit(1)
        except Exception as e:
            print(e)
            sys.exit(1)

    @staticmethod
    def WriteResultToFile(pair_out_filename, res):
        with open(pair_out_filename, 'w+') as f:
            for line in res:
                f.write(line+'\n')

def PrintReport(syst, sapt, scsig, polarisable):
    report = "Epsilon LJ parameters were scaled by "
    if sapt:
        report += "k_sapt"
    else:
        report += "k_pred"
    report += " parameter"

    if polarisable is not None:
        tmp = ', '.join(polarisable)
        report += ". Fragments %s were already polarisable " % tmp

    report +=". Changes are marked with '~'.\n"

    report += "Sigma LJ parameters "
    if scsig is None:
        report += "were not scaled.\n"
    else: 
        report +="were scaled by %5.3f value." %  _Const.sigma_k
        
        if all([x.scale_sig for x in syst.fragments]):
            report += " All fragments were scaled."
        else:
            sig_list = [f.name for f in syst.fragments if f.scale_sig]
            tmp = ', '.join(sig_list)
            report += " Only %s fragments were scaled." % tmp
        
        report +=" Changes are marked with '*'.\n"

    report += '------------------------------------------\n'

    report += 'Fragment i   Fragment j   k_sapt'
    if not sapt:
        report += '    k_pred'

    for frp in syst.fragmentpairs:
        report += '\n%-10s   %-10s' % (frp.fr1.name, frp.fr2.name)
        if sapt:
            report += '  %6.2f' % frp.k_sapt
        elif frp.k_sapt is None:
            report += '   %6s    %6.2f' % ('-',frp.k_pred)
        else:
            report += '   %6.2f    %6.2f' % (frp.k_sapt,frp.k_pred)

    report += '\n------------------------------------------'
    print(report)

def main():
    parser = argparse.ArgumentParser(description = usage, formatter_class = argparse.RawTextHelpFormatter)
    parser.add_argument('-f', '--ff_filename', type=str, default = 'fragment.ff', help = 'fragment force field (default: fragment.ff)')
    parser.add_argument('-a', '--alpha_filename', type=str, default = 'alpha.ff', help = 'polarisability values file (default: alpha.ff)')
    parser.add_argument('-i', '--input_filename', type=str, default = 'fragment.inp', help = 'fragment input file with atomic indices (default: fragment.inp)')
    parser.add_argument('-ip', '--pair_in_filename', type=str, default = 'pair-p.lmp', help = 'pair style input file (default: pair-p.lmp)')
    parser.add_argument('-op', '--pair_out_filename', type=str, default = 'pair-p-sc.lmp', help = 'pair style output file (default: pair-p-sc.lmp)')
    parser.add_argument('-q', '--sapt', action = 'store_true', help = 'use sapt calculated k values, default: use predicted k values')
    parser.add_argument('-s', '--scsig', nargs='*', type = str,  help = 'scale sigma if specified; default value: 0.985; \n\
    -s                       - scale all fragments\' sigma by 0.985 \n\
    -s value                 - scale all fragments\' sigma by user-defined value \n\
    -s name1 name2 ...       - scale the specified fragments\' sigma by 0.985 \n\
    -s value name1 name2 ... - scale the specified fragments\' sigma by user-defined value')
    parser.add_argument('-p', '--polarisable', nargs='+', type = str, help = 'already polarisable monomers')

    args = parser.parse_args()

    ff = Forcefield(args.ff_filename)
    syst = System(args.input_filename)
    syst.GetFragments(ff,args.polarisable)

    if args.scsig is not None:
        syst.ParseScaleSigma(args.scsig)

    syst.GetFragmentPairs(ff)

    if (not args.sapt):
        pol = Polarisation(args.alpha_filename)
        syst.GetFragAtoms(pol)
        syst.GetKPred()
    
    res = ScaleLJ.Scale(args.pair_in_filename,ff,syst,args.sapt,args.scsig)
    ScaleLJ.WriteResultToFile(args.pair_out_filename, res)

    PrintReport(syst,args.sapt,args.scsig,args.polarisable)
    
if __name__ == '__main__':
    main()