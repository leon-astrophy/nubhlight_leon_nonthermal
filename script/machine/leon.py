################################################################################
#                                                                              #
#  MACHINE-SPECIFIC FUNCTIONS                                                  #
#                                                                              #
#    OPTIONS:                                                                  #
#      COMPILER   : PATH TO COMPILER EXECUTABLE                                #
#      GSL_DIR    : PATH TO GSL INSTALLATION                                   #
#      MPI_DIR    : PATH TO MPI INSTALLATION                                   #
#      HDF5_DIR   : PATH TO HDF5 INSTALLATION                                  #
#      EXECUTABLE : BINARY WRAPPER USED TO LAUNCH BHLIGHT                      #
#                                                                              #
#    MPI_DIR AND HDF5_DIR ARE NOT REQUIRED IF COMPILER HANDLES HEADERS AND     #
#    LIBRARIES FOR THESE DEPENDENCIES                                          #
#                                                                              #
################################################################################

# import #
import os

# compilation flag, copy from iharm3d #
flags_base = '-std=gnu99 -O3 -march=native -mtune=native -flto -funroll-loops -ftree-vectorize -pipe -fopenmp'

# host machine name #
def matches_host():
  host = os.uname()[1]
  return host == 'gpus'

# get host attribute #
def get_options():
  host = {}
  host['NAME']           = os.uname()[1]
  host['COMPILER']       = 'h5pcc'
  host['COMPILER_FLAGS'] = flags_base + ' ' + ' '
  host['DEBUG_FLAGS']    = flags_base + ' ' + ' '
  host['GSL_DIR']        = '/usr/include/gsl'

  # return #
  return host