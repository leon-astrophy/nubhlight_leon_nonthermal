################################################################################
#                                                                              #
#  OPTICALLY THIN BREMSSTRAHLUNG COOLING                                       #
#                                                                              #
################################################################################

# import #
import sys; sys.path.append('../../script/'); 
sys.dont_write_bytecode = True; import bhlight as bhl; del sys

################################################################################

# problem name #
PROB = 'brem'

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
bhl.config.set_cparm('X1L_GAS_BOUND', 'BC_OUTFLOW')
bhl.config.set_cparm('X1R_GAS_BOUND', 'BC_OUTFLOW')
bhl.config.set_cparm('X2L_GAS_BOUND', 'BC_OUTFLOW')
bhl.config.set_cparm('X2R_GAS_BOUND', 'BC_OUTFLOW')
bhl.config.set_cparm('X3L_GAS_BOUND', 'BC_OUTFLOW')
bhl.config.set_cparm('X3R_GAS_BOUND', 'BC_OUTFLOW')
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
bhl.config.set_cparm('EMISSION', True)
bhl.config.set_cparm('ABSORPTION', False)
bhl.config.set_cparm('SCATTERING', False)
bhl.config.set_cparm('SYNCHROTRON', False)
bhl.config.set_cparm('BREMSSTRAHLUNG', True)
bhl.config.set_cparm('GRAYABSORPTION', False)
bhl.config.set_cparm('EXPTAU_WEIGHTS', True)
bhl.config.set_cparm('ESTIMATE_THETAE', False)
bhl.config.set_cparm('X1L_RAD_BOUND', 'BC_ESCAPE')
bhl.config.set_cparm('X1R_RAD_BOUND', 'BC_ESCAPE')
bhl.config.set_cparm('X2L_RAD_BOUND', 'BC_ESCAPE')
bhl.config.set_cparm('X2R_RAD_BOUND', 'BC_ESCAPE')
bhl.config.set_cparm('X3L_RAD_BOUND', 'BC_ESCAPE')
bhl.config.set_cparm('X3R_RAD_BOUND', 'BC_ESCAPE')

################################################################################
### RUNTIME PARAMETERS ###

# simulation time #
bhl.config.set_rparm('tf', 'double', default = 256.e-0)
bhl.config.set_rparm('dt', 'double', default = 256.e-6)

# physical unit #
bhl.config.set_rparm('L_unit', 'double', default = 1.17106428906e16)
bhl.config.set_rparm('M_unit', 'double', default = 1.57378000000e32)

# output frequency #
bhl.config.set_rparm('DTd', 'double', default = 5.0e-1)
bhl.config.set_rparm('DTl', 'double', default = 10.e-0)
bhl.config.set_rparm('DTr', 'integer', default = 5000000.)
bhl.config.set_rparm('DNr', 'integer', default = 5000000.)

# radiation #
bhl.config.set_rparm('tune_scatt', 'double', default = 0.00)
bhl.config.set_rparm('tune_emiss', 'double', default = 100.)
bhl.config.set_rparm('t0_tune_emiss', 'double', -1)
bhl.config.set_rparm('t0_tune_scatt', 'double', -1)

# extra parameters, initial temperature #
bhl.config.set_rparm('T0', 'double', default = 1.e8)

################################################################################
### CONFIGURE AND COMPILE  ###

bhl.build(PROB)

