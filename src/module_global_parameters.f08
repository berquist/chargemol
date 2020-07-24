 MODULE module_global_parameters
 !===================================================================================
 ! Global parameters.
 ! The default values should give good results for all chemical systems, so
 ! don't change them unless you have a compelling reason.
 ! Copyright (c) 2014, 2015, 2016, 2017 by Thomas A. Manz and Nidia Gabaldon Limas. Rights reserved.
 !===================================================================================
 USE module_precision
 
 IMPLICIT NONE
 
 INTEGER, PARAMETER :: nshells = 100 !Number of radial integration shells
 INTEGER, PARAMETER :: num_lookup_points = 10000
 INTEGER, PARAMETER :: read_buffer_size = 21621600 !This number is evenly divisible by all whole numbers from 1 to 16 (inclusive).
!$ INTEGER, PARAMETER :: collapsed_chunk_size = 2500 
 REAL(kind=dp), PARAMETER :: Xi_threshold = 0.001_dp
 REAL(kind=dp), PARAMETER :: pi = 3.14159265358979323846_dp
 REAL(kind=dp), PARAMETER :: Xi_zero_tolerance = 0.000001_dp 
 REAL(kind=dp), PARAMETER :: cutoff_radius = 500.0_dp !picometers
 REAL(kind=dp), PARAMETER :: bohrperangstrom = 1.889725989_dp
 REAL(kind=dp), PARAMETER :: zero_tolerance = 1.0e-10_dp 
 REAL(kind=dp), PARAMETER :: integration_tolerance = 0.10_dp
 REAL(kind=dp), PARAMETER :: integration_tolerance_percent = 0.10_dp
 REAL(kind=dp), PARAMETER :: maxpixelvolume = 0.0157_dp
 REAL(kind=dp), PARAMETER :: pixel_integration_tolerance = 0.03_dp
 REAL(kind=dp), PARAMETER :: charge_convergence_tolerance = 1.0e-5_dp
 REAL(kind=dp), PARAMETER :: spin_convergence_tolerance = 5.0e-5_dp
 REAL(kind=dp), PARAMETER :: distance_tolerance = 1.0e-6_dp
 REAL(kind=dp), PARAMETER :: rmin_cloud_penetration=200.0_dp !picometers
 REAL(kind=dp), PARAMETER :: spin_ref_fraction = 1.0_dp/2.0_dp
 REAL(kind=dp), PARAMETER :: wA_renormalization_max=10.0_dp
 REAL(kind=dp), PARAMETER :: scalefactor = 100.0_dp*nshells/(cutoff_radius*bohrperangstrom) !10.58354497764173_dp !
 REAL(kind=dp), PARAMETER :: localizing_power = 4.0_dp
 REAL(kind=dp), PARAMETER :: max_core_electrons_per_pixel = 30.0_dp
 !Parameters used only for computing the effective bond orders (BOs)
 REAL(kind=dp), PARAMETER :: BO_print_cutoff = 0.001_dp !BOs less than this value are not printed
 !Charge computation parameters
 REAL(kind=dp), PARAMETER :: minimum_buried_tail_exponent = 1.75_dp
 REAL(kind=dp), PARAMETER :: maximum_buried_tail_exponent = 2.5_dp
 !These are the standard atomic weights of the elements as taken from the NIST.GOV website
 !and published by T.B. Coplen1 in Atomic Weights of the Elements 1999
 !they include changes reported from the 2001 review in Chem. Int., 23, 179 (2001))
 !these are averaged over natural isotope abundances
 REAL(kind=dp), PARAMETER,DIMENSION(109) :: atomic_weight=[1.00794_dp, 4.002602_dp, 6.941_dp,9.012182_dp,10.811_dp,12.0107_dp,&
 14.0067_dp,15.9994_dp,18.9984032_dp,20.1797_dp,22.989770_dp,24.3050_dp,26.981538_dp,28.0855_dp,30.973761_dp,32.065_dp,35.453_dp,&
 39.948_dp,39.0983_dp,40.078_dp,44.955910_dp,47.867_dp,50.9415_dp,51.9961_dp,54.938049_dp,55.845_dp,58.933200_dp,58.6934_dp,&
 63.546_dp,65.409_dp,69.723_dp,72.64_dp,74.92160_dp,78.96_dp,79.904_dp,83.798_dp,85.4678_dp,87.62_dp,88.90585_dp,91.224_dp,&
 92.90638_dp,95.94_dp,98.0_dp,101.07_dp,102.90550_dp,106.42_dp,107.8682_dp,112.411_dp,114.818_dp,118.710_dp,121.760_dp,&
 127.60_dp,126.90447_dp,131.293_dp,132.90545_dp,137.327_dp,138.9055_dp,140.116_dp,140.90765_dp,144.24_dp,145.0_dp,150.36_dp,&
 151.964_dp,157.25_dp,158.92534_dp,162.500_dp,164.93032_dp,167.259_dp,168.93421_dp,173.04_dp,174.967_dp,178.49_dp,180.9479_dp,&
 183.84_dp,186.207_dp,190.23_dp,192.217_dp,195.078_dp,196.96655_dp,200.59_dp,204.3833_dp,207.2_dp,208.98038_dp,209.0_dp,210.0_dp,&
 222.0_dp,223.0_dp,226.0_dp,227.0_dp,232.0381_dp,231.03588_dp,238.02891_dp,237.0_dp,244.0_dp,243.0_dp,247.0_dp,247.0_dp,251.0_dp,&
 252.0_dp,257.0_dp,258.0_dp,259.0_dp,262.0_dp,261.0_dp,262.0_dp,266.0_dp,264.0_dp,277.0_dp,268.0_dp]
 ! Default direction of SAXIS in VASP
 REAL(kind=dp), DIMENSION(3) :: SAXIS = (/1.0e-12_dp, 0.0_dp, 1.0_dp/)
 ! Parameters used only for reading gaussian basis set wfx files
 REAL(kind=dp), PARAMETER :: preferred_grid_spacing = 0.14_dp
 REAL(kind=dp), PARAMETER :: gaussian_overlap_tolerance = 1.0e-12_dp
 REAL(kind=dp), PARAMETER :: analytic_alpha_cutoff = 5.0_dp !primitive pairs with alpha sum greater than this are treated analytically
 REAL(kind=dp), PARAMETER :: periodic_cutoff_length = 28.0_dp !This is in bohr
 REAL(kind=dp), PARAMETER :: coarser_grid_alpha_cutoff = 0.4_dp !primitive pairs with alpha sum less than this area integrated over coarser grids
 
 ! Some default input files for different types of jobs
 ! Default name for xsf file
 CHARACTER(19), PARAMETER :: xsf_inputfile = 'valence_density.xsf'
 !Version information
 CHARACTER(41), PARAMETER :: version = 'Chargemol version 3.5 September 26, 2017.'
 CHARACTER(2), PARAMETER :: density_set_prefix = 'c2'
 
 LOGICAL :: run_interactive=.FALSE.
 
 END MODULE module_global_parameters