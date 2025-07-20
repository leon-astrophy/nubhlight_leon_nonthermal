################################################################################
#                                                                              #
#  COMPTONIZATION                                                              #
#                                                                              #
################################################################################

# import #
import sys; sys.path.append('../../script/'); 
sys.dont_write_bytecode = True; import bhlight as bhl; del sys

################################################################################

# problem name #
PROB = 'comptonization'

################################################################################
### COMPILE TIME PARAMETERS ###

# SPATIAL RESOLUTION AND MPI DECOMPOSITION
bhl.config.set_cparm('N1TOT', 1)
bhl.config.set_cparm('N2TOT', 1)
bhl.config.set_cparm('N3TOT', 1)
bhl.config.set_cparm('N1CPU', 1)
bhl.config.set_cparm('N2CPU', 1)
bhl.config.set_cparm('N3CPU', 1)

# OPENMP PARALLELIZATION
bhl.config.set_cparm('OPENMP', True)

# COORDINATES
bhl.config.set_cparm('METRIC', 'MINKOWSKI')

# FLUID
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

# RADIATION
bhl.config.set_cparm('NTH', 8)
bhl.config.set_cparm('NPHI', 8)
bhl.config.set_cparm('NU_BINS', 200)
bhl.config.set_cparm('RADIATION', True)
bhl.config.set_cparm('EMISSION', False)
bhl.config.set_cparm('ABSORPTION', False)
bhl.config.set_cparm('SCATTERING', True)
bhl.config.set_cparm('SYNCHROTRON', False)
bhl.config.set_cparm('BREMSSTRAHLUNG', False)
bhl.config.set_cparm('GRAYABSORPTION', False)
bhl.config.set_cparm('EXPTAU_WEIGHTS', True)
bhl.config.set_cparm('ESTIMATE_THETAE', False)
bhl.config.set_cparm('X1L_RAD_BOUND', 'BC_PERIODIC')
bhl.config.set_cparm('X1R_RAD_BOUND', 'BC_PERIODIC')
bhl.config.set_cparm('X2L_RAD_BOUND', 'BC_PERIODIC')
bhl.config.set_cparm('X2R_RAD_BOUND', 'BC_PERIODIC')
bhl.config.set_cparm('X3L_RAD_BOUND', 'BC_PERIODIC')
bhl.config.set_cparm('X3R_RAD_BOUND', 'BC_PERIODIC')

################################################################################
### RUNTIME PARAMETERS ###

# simulation time #
bhl.config.set_rparm('tf', 'double', default = 10000.e-0)
bhl.config.set_rparm('dt', 'double', default = 10000.e-6)

# physical unit #
bhl.config.set_rparm('L_unit', 'double', default = 1.503200e6)
bhl.config.set_rparm('M_unit', 'double', default = 1.42034e12)

# output frequency #
bhl.config.set_rparm('DTd', 'double', default = 10.)
bhl.config.set_rparm('DTl', 'double', default = 500.)
bhl.config.set_rparm('DTr', 'integer', default = 5000000.)
bhl.config.set_rparm('DNr', 'integer', default = 5000000.)

# radiation #
bhl.config.set_rparm('tune_scatt', 'double', default = 0.)
bhl.config.set_rparm('tune_emiss', 'double', default = 2.)
bhl.config.set_rparm('t0_tune_emiss', 'double', -1)
bhl.config.set_rparm('t0_tune_scatt', 'double', -1)

################################################################################
### CONFIGURE AND COMPILE  ###

bhl.build(PROB)

