################################################################################
#                                                                              #
#  Nonthermal problem in 1D                                                    #
#                                                                              #
################################################################################

# import #
import sys; sys.path.append('../../script/')
sys.dont_write_bytecode = True; import bhlight as bhl;

################################################################################

# problem name #
PROB = 'nonthermal_synchrotron'

################################################################################
### COMPILE TIME PARAMETERS ###

# SPATIAL RESOLUTION AND MPI DECOMPOSITION #
bhl.config.set_cparm('N1TOT', 1)
bhl.config.set_cparm('N2TOT', 1)
bhl.config.set_cparm('N3TOT', 1)
bhl.config.set_cparm('N1CPU', 1)
bhl.config.set_cparm('N2CPU', 1)
bhl.config.set_cparm('N3CPU', 1)

# OPENMP PARALLELIZATION #
bhl.config.set_cparm('OPENMP', True)

# COORDINATES #
bhl.config.set_cparm('METRIC', 'MINKOWSKI')

# ELECTRONS #
bhl.config.set_cparm('ELECTRONS', True)
bhl.config.set_cparm('BETA_HEAT', False)
bhl.config.set_cparm('COULOMB', False)
bhl.config.set_cparm('SUPPRESS_HIGHB_HEAT', False)

# Do log dumping #
bhl.config.set_cparm('LOGDUMPING', True)

# NONTHERMAL #
bhl.config.set_cparm('NONTHERMAL', True)
bhl.config.set_cparm('LIMIT_NTE', False)
bhl.config.set_cparm('COMPTON_NTE', False)
bhl.config.set_cparm('COULOMB_NTE', False)
bhl.config.set_cparm('SYNCHROTRON_NTE', True)
bhl.config.set_cparm('BREMSSTRAHLUNG_NTE', False)
bhl.config.set_cparm('ADIABTIC_SCALING', False)
bhl.config.set_cparm('PLAW', 3.5)

### SYNCHROTRON/ADIABATIC PARAMS ###
bhl.config.set_cparm('ART_SYNCH', 200)
bhl.config.set_cparm('const_inj', 1.0e0)
bhl.config.set_cparm('ART_ADIAB', 0.0e0)
bhl.config.set_cparm('SKIP_ADIAB', True)
bhl.config.set_cparm('SKIP_VISCOUS', False)
bhl.config.set_cparm('SKIP_COOLING', False)
bhl.config.set_cparm('CONST_INJECTION', True)

# FLUID #
bhl.config.set_cparm('RECONSTRUCTION', 'WENO')
bhl.config.set_cparm('X1L_GAS_BOUND', 'BC_PERIODIC')
bhl.config.set_cparm('X1R_GAS_BOUND', 'BC_PERIODIC')
bhl.config.set_cparm('X2L_GAS_BOUND', 'BC_PERIODIC')
bhl.config.set_cparm('X2R_GAS_BOUND', 'BC_PERIODIC')
bhl.config.set_cparm('X3L_GAS_BOUND', 'BC_PERIODIC')
bhl.config.set_cparm('X3R_GAS_BOUND', 'BC_PERIODIC')
bhl.config.set_cparm('X1L_INFLOW', False)
bhl.config.set_cparm('X1R_INFLOW', False)
bhl.config.set_cparm('X2L_INFLOW', False)
bhl.config.set_cparm('X2R_INFLOW', False)
bhl.config.set_cparm('X3L_INFLOW', False)
bhl.config.set_cparm('X3R_INFLOW', False)

################################################################################
### RUNTIME PARAMETERS ###

# gas parameters #
bhl.config.set_rparm('gam', 'double', default = 1.4)

# simulation time #
bhl.config.set_rparm('tf', 'double', default = 10000)
bhl.config.set_rparm('dt', 'double', default = 0.0001)

# simulation unit #
bhl.config.set_rparm('T_unit', 'double', default = 1.0e0)
bhl.config.set_rparm('L_unit', 'double', default = 1.0e0)
bhl.config.set_rparm('M_unit', 'double', default = 1.0e0)

# output frequency #
bhl.config.set_rparm('DTd', 'double', default = 0.001)
bhl.config.set_rparm('DTl', 'double', default = 1000000)
bhl.config.set_rparm('DTr', 'double', default = 1000000)
bhl.config.set_rparm('DNr', 'integer', default = 1000000)

################################################################################
### CONFIGURE AND COMPILE  ###

bhl.build(PROB)
