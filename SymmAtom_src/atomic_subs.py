import numpy as np
from atomic_constants import * 
from scipy.special import comb, perm
from timeit import default_timer as timer
from datetime import timedelta, datetime
from time import gmtime, strftime

def show_header():
    print(11*' '+'*'+54*'-'+'*')
    print(11*' '+'|'+54*' '+'|')
    print(11*' '+'|'+23*' '+'SymmAtom'+23*' '+'|')
    print(11*' '+'|'+54*' '+'|')
    print(11*' '+'*'+54*'-'+'*')
    print(11*' '+'|'+54*' '+'|')
    wel = 'welcome to Symmetry labeled atom subroutine'
    print('{}{:<3}{:^50}{:>3}'.format(11*' ','|',wel,'|'))
    gitstr="git@github.com:alex2shiyu/atom2020.git"
    print('{}{:<3}{:^50}{:>3}'.format(11*' ','|',gitstr,'|'))
    print(11*' '+'|'+54*' '+'|')
    print(11*' '+'|'+54*' '+'|')
    print('{}{:<3}{:<50}{:>3}'.format(11*' ','|','Authors: ','|'))
    print('{}{:<3}{:<50}{:>3}'.format(11*' ','|','       : Shiyu Peng (sypeng@iphy.ac.cn)','|'))
    print('{}{:<3}{:<50}{:>3}'.format(11*' ','|','       : Xi Dai','|'))
    print(11*' '+'|'+54*' '+'|')
    print('{}{:<3}{:<50}{:>3}'.format(11*' ','|','This program is distributed in the hope that it','|'))
    print('{}{:<3}{:<50}{:>3}'.format(11*' ','|','will be useful, but WITHOUT ANY WARRANTY; without','|'))
    print('{}{:<3}{:<50}{:>3}'.format(11*' ','|','even the implied warranty of MERCHANTABILITY or','|'))
    print('{}{:<3}{:<50}{:>3}'.format(11*' ','|','FITNESS FOR A PARTICULAR PURPOSE. See the GNU','|'))
    print('{}{:<3}{:<50}{:>3}'.format(11*' ','|','General Public License for more details.','|'))
    print(11*' '+'*'+54*'-'+'*')
    print('{}{:<3}{:^50}{:>3}'.format(11*' ','|','Release: 1.0        4th June 2020','|'))
    print(11*' '+'*'+54*'-'+'*')
    print('\n\n')


def show_subheader(header):
    print('\n\n')
    print(11*' '+'*'+54*'='+'*')
    print('{}{:<3}{:^50}{:>3}'.format(11*' ','|','>>>'+header+'<<<','|'))
    print(11*' '+'*'+54*'='+'*')

def show_subsubheader(header):
    print('\n')
    print(18*' ',40*'-')
    print('{}{}{:^38}{}'.format(19*' ','*',header,'*'))
    print(18*' ',40*'-')
    print('\n')

def show_sub3header(header):
    print('\n')
    print(30*'-')
    print('{}{:<40}'.format('|>      ',header))
    print(30*'-')

def show_sub4header(header):
    print('\n')
    print(20*'-')
    print('{}{:<40}'.format('|>      ',header))
    print(20*'-')

def show_error(header):
    print('\n')
    print('*'+24*'='+'*')
    print('{:<3}{:^20}{:>3}'.format('|','Error in'+header,'|'))
    print('*'+24*'='+'*')

def show_success():
    print('\n')
#   print(20*'*')
    print('**   success ^_^  **')
#   print(20*'*')


def show_over():
    print('\n\n')
    print(' *'+75*'*'+'*')
    print(' *'+31*' ','All Done !',32*' '+'*')
    print(' *'+75*'*'+'*')

def atomic_state(norbs):
    from scipy.special import comb
    nstat = np.zeros(norbs+1,dtype=np.int32)
    for i in range(norbs+1):
        nstat[i] = comb(norbs,i)
    return nstat

def decide_unitary_len(nocc,unitary_mat):
    '''
    aim : to decide the unitary_len needed by fortran warperred function gw_make_newui
    with the precision cutoff setting
    intput : 
            nocc : [integer] to decide the subspace
            unitary : [numpy.ndarray] the unitary matrix in single particle basis
    '''
    norb = unitary_mat.shape[0]
    norb_eff = int(0)
    norb_eff_set = []
    for iorb in range(norb):
        norb_eff_set[iorb] = sum(np.abs(unitary_mat[:,iorb]) > unitary_manybody_prec)
    norb_eff = int(np.max(norb_eff_set))

    unitary_len = perm(norb_eff,nocc) + int(100)
    return unitary_len

class atom_timer():
    '''
    aim : counting the elapsed time of diferent block and functions
    '''
    def __init__(self):
        self.start = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
        self.atom  = 0.0
        self.cumat = 0.0
        self.basis = 0.0
        self.ham  = 0.0
        self.dpg   = 0.0
        self.TranOrb = 0.0
        self.Nsubstot = 0.0
        self.dump = 0.0
        self.collect = 0.0
        self.Nsubs = []
        self.before = 0.0
        self.end   = 0.0
    
    def reset(self):
        self.before = timer()

    def count(self):
        time_now = timer()
        # here we assume the elapsed time less than 1 day
#       return str(timedelta(seconds=(time_now - self.before)))[0:7]
        return str(timedelta(seconds=(time_now - self.before)))[0:14]

    def timeover(self):
        self.end = datetime.now().strftime('%Y-%m-%d %H:%M:%S')

    
    def show(self):
        print('\n\n')
        print(' *'+75*'='+'*')
        print(' |'+29*' ','TIME INFORMATION',28*' '+'|')
        print(' *'+75*'='+'*')
        print(' |'+'  Tag   ',57*' ','Time(s)','|')
        print(' |'+75*'-'+'|')
        print('{:>2}{:<20}{:>31}{:>23}{:>2}'.format('|','  start',':',self.start,'|')) 
        print('{:>2}{:<20}{:>31}{:>23}{:>2}'.format('|','  atom',':',self.atom,'|')) 
        print('{:>2}{:<20}{:>31}{:>23}{:>2}'.format('|','  cumat',':',self.cumat,'|')) 
        print('{:>2}{:<20}{:>31}{:>23}{:>2}'.format('|','  basis',':',self.basis,'|')) 
        print('{:>2}{:<20}{:>31}{:>23}{:>2}'.format('|','  hamiltonain',':',self.ham,'|')) 
        print('{:>2}{:<20}{:>31}{:>23}{:>2}'.format('|','  double PG',':',self.dpg,'|')) 
        print('{:>2}{:<20}{:>31}{:>23}{:>2}'.format('|','  TranOrb',':',self.TranOrb,'|')) 
        print('{:>2}{:<20}{:>31}{:>23}{:>2}'.format('|','  Total of Nsubs',':',self.Nsubstot,'|')) 
        for itime in self.Nsubs :
            itime.show()
        print('{:>2}{:<20}{:>31}{:>23}{:>2}'.format('|','  collect',':',self.collect,'|')) 
        print('{:>2}{:<20}{:>31}{:>23}{:>2}'.format('|','  dump',':',self.dump,'|')) 
        print('{:>2}{:<20}{:>31}{:>23}{:>2}'.format('|','  end',':',self.end,'|')) 
        print(' *'+75*'='+'*')


class atom_timer_Nsubs(atom_timer):
    '''
    aim : counting time characters of different sub Fock space of N
    '''
    def __init__(self):
        super().__init__()
        self.operators = 0
        self.check_ham = 0
        self.check_pro = 0
        self.reduct    = 0
        self.tran2natural = 0
        self.vpm = 0
        self.decompose = 0
        self.check_irrepbasis = 0
        self.noc       = 0
        self.ham_trandiag = 0
        self.check_eigenwave_irrep = 0
        
    def show(self):
        print('{:>2}{}{:<20}{:>25}{:>23}{:>2}'.format('|',6*' ','N = '+str(self.noc),' ',' ','|')) 
        print('{:>2}{}{:<25}{:>15}{:>23}{:>2}'.format('|',11*' ','basis',':',self.basis,'|')) 
        print('{:>2}{}{:<25}{:>15}{:>23}{:>2}'.format('|',11*' ','operators',':',self.operators,'|')) 
        print('{:>2}{}{:<25}{:>15}{:>23}{:>2}'.format('|',11*' ','check ham',':',self.check_ham,'|')) 
        print('{:>2}{}{:<25}{:>15}{:>23}{:>2}'.format('|',11*' ','reduction',':',self.reduct,'|')) 
        print('{:>2}{}{:<25}{:>15}{:>23}{:>2}'.format('|',11*' ','check projections',':',self.check_pro,'|')) 
        print('{:>2}{}{:<25}{:>15}{:>23}{:>2}'.format('|',11*' ','check PG',':',self.check_irrepbasis,'|')) 
        print('{:>2}{}{:<25}{:>15}{:>23}{:>2}'.format('|',11*' ','trans, diag ham',':',self.ham_trandiag,'|')) 
        print('{:>2}{}{:<25}{:>15}{:>23}{:>2}'.format('|',11*' ','check eigenwave(symm)',':',self.check_eigenwave_irrep,'|')) 
        print('{:>2}{}{:<25}{:>15}{:>23}{:>2}'.format('|',11*' ','decomposition',':',self.decompose,'|')) 
        print('{:>2}{}{:<25}{:>15}{:>23}{:>2}'.format('|',11*' ','trans to natural',':',self.tran2natural,'|')) 
        print('{:>2}{}{:<25}{:>15}{:>23}{:>2}'.format('|',11*' ','make vpms',':',self.vpm,'|')) 
