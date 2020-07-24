 MODULE module_common_variable_declarations
 !=====================================================================================
 ! Variables used by all files
 ! Copyright (c) 2014, 2015, 2016, 2017 by Thomas A. Manz and Nidia Gabaldon Limas. Rights reserved.
 !=====================================================================================

 USE module_precision
 !$ USE omp_lib
 
 IMPLICIT NONE

 !-----------------------------------------------------------------------------------
 !Public variables
 !-----------------------------------------------------------------------------------

 CHARACTER(3) :: input_type
 CHARACTER(5)  :: zone
 CHARACTER(8)  :: date,charge_type
 CHARACTER(10) :: time
 CHARACTER(11) :: density_format
 CHARACTER(200) :: atomic_densities_directory,output_directory_path
 CHARACTER(200):: basis_set_type,long_input_filename,job_directory
 CHARACTER(:), ALLOCATABLE :: output_filename
 INTEGER,DIMENSION(8) :: values
 INTEGER :: i,j,flag,num_periodic_directions,base_length,totnumA,totnumB,totnumC,delta_na,delta_nb,&
 delta_nc,ncenters,nprimitives,n_edf_primitives,norbitals,included_electrons,input_FID,output_FID,natoms,&
 na,nb,nc,ka,kb,kc,iter,spin_iter,potcar_FID,iostat_value,total_density_FID,valence_FID,nsingleblocks,start_index,stop_index,&
 shell_index,ka_plus,kb_plus,kc_plus,number_of_threads,chunk_size,num_corrected_pixels,num_sym_unique_bond_pairs,&
 num_nonzero_components,active_pair,BO_FID,current_pair,job_control_FID,&
 available_reference_ion_range(109,2),radial_moment_FID,num_core(109),periA,periB,periC,buffer
 INTEGER, ALLOCATABLE, DIMENSION (:) :: atomic_number,Apoints,Bpoints,Cpoints,missing_core,&
 nprimitives_per_center,lower_na,upper_na,lower_nb,upper_nb,lower_nc,upper_nc
 INTEGER,ALLOCATABLE, DIMENSION (:,:) :: center_nabc,sum_points,single_block_ranges,kApoints,kBpoints,kCpoints
 INTEGER,ALLOCATABLE,DIMENSION(:,:,:) :: dominant_atom_points
 REAL(kind=dp) :: netcharge,int_tolerance,boundary(3,3),origin(3),pixelvolume,vector1(3),vector2(3),vector3(3),&
 max_change,old_change,old_old_change,old_old_old_change,ncore,nvalence,normalization,center_of_mass(3),checkme,&
 tot_spin_moment_vector(3),tot_spin_moment,inv_boundary(3,3),seconds,max_correction_size,MAD1,&
 distance_from_nucleus,partial_charge_dipole_x,&
 partial_charge_dipole_y,partial_charge_dipole_z,partial_charge_dipole_magnitude,partial_charge_quadrupole_xy,&
 partial_charge_quadrupole_xz,partial_charge_quadrupole_yz,partial_charge_quadrupole_x2minusy2,sum_negative_density,&
 partial_charge_quadrupole_3z2minusr2,partial_charge_quadrupole_eigenvals(3),dipole_x,dipole_y,dipole_z,dipole_magnitude,&
 quadrupole_xy,quadrupole_xz,quadrupole_yz,quadrupole_x2minusy2,quadrupole_3z2minusr2,quadrupole_eigenvals(3),distance_scale,&
 full_quadrupole_xy,full_quadrupole_xz,full_quadrupole_yz,full_quadrupole_x2minusy2,valence_electrons,reference_weighting,&
 full_quadrupole_3z2minusr2,full_quadrupole_eigenvals(3),parameters(4,4),density_scale
 REAL(kind=dp), ALLOCATABLE, DIMENSION(:) :: valence_population,core_population,change,core_electrons,&
 fitted_tail_slope,fitted_tail_intercept,fitted_tail_Rsquared, spin_population,net_atomic_charge,&
 radial_shell_volume,integrated_nA,summed_contact_exchange,max_density,K_factor,&
 reference_ion_charge,Hirshfeld_net_atomic_charge,localized_net_atomic_charge,localized_population,&
 old_net_atomic_charge,old_reference_ion_charge,atomic_polarizability_upper_bound,CCSD_free_atom_polarizability_to_volume_ratio
 REAL(kind=dp),ALLOCATABLE,DIMENSION(:) :: effective_nuclear_charge,atomic_number2,effective_nuclear_charge2,atomic_summed_BO
 REAL(kind=dp), ALLOCATABLE, DIMENSION(:,:) :: basis_set_centers,primitive_exponents, alpha_beta_occupation,orbital_coefficients,&
 edf_primitives,center_shift,occupancy_correction,periodic_vectors,coords,spherical_average_density,&
 spherical_avg_ref_pseudodensity,avg_PtoWref,neutral_density,combined_oxidation_density,partial_density,&
 conditioned_oxidation_density,partial_core_density,final_result,spherical_avg_core_density,sum_core_density,sum_density,&
 spin_population_vector,dominant_atom_weight,bond_pair_matrix,corrected_spherical_avg_density,nA_outside,sum_conditioned_density,&
 spherical_average_atomic_spin,tau,coords2,raw_matrix_cube,min_density,shell_analytic_radial_moments
 REAL(kind=dp), ALLOCATABLE, DIMENSION (:,:,:) :: valence_density,spin_density,core_density,total_pseudodensity,&
 core_pseudodensity,corrected_total_density,local_spherical_avg_atomic_exchange_vectors,localized_pseudodensity,&
 dot_product_total_spherical_avg_atomic_exchange_vectors,total_density,spherical_average_atomic_spin_vector
 REAL(kind=dp), ALLOCATABLE, DIMENSION (:,:,:,:) :: spin_density_vector,total_local_spherical_avg_atomic_exchange_vectors
 LOGICAL :: spin_available,core_available, non_collinear,periodicA, periodicB, periodicC,edf_available,run_parallel,&
 atomic_reference_polarizabilities_available,compute_BOs,print_atomic_densities
 
 END MODULE module_common_variable_declarations