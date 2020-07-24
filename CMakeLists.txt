cmake_minimum_required(VERSION 3.1)
project(chargemol Fortran)

option(ENABLE_OPENMP "Enable OpenMP parallelization" ON)
option(ENABLE_STATIC_LINK "Force static linkage" OFF)

set(SRC
    src/module_precision.f08
    src/module_global_parameters.f08
    src/module_common_variable_declarations.f08
    src/module_quote.f08
    src/module_string_utilities.f08
    src/module_matrix_operations.f08
    src/module_reshaping_functions.f08
    src/module_CM5_parameters.f08
    src/module_calculate_atomic_polarizabilities_upper_bound.f08
    src/module_read_job_control.f08
    src/module_read_wfx.f08
    src/module_gaussian_functions.f08
    src/module_compute_CM5.f08
    src/module_check_atomic_reference_polarizabilities_available.f08
    src/module_charge_center_positions_and_parallelpiped.f08
    src/module_generate_kApoints_kBpoints_kCpoints.f08
    src/module_check_grid_spacing.f08
    src/module_add_missing_core_density.f08
    src/module_oxidation_density.f08
    src/module_initialize_atomic_densities.f08
    src/module_align_periodic_atoms.f08
    src/module_read_spin_density_cube_files.f08
    src/module_format_valence_and_total_cube_densities.f08
    src/module_format_valence_cube_density.f08
    src/module_format_total_cube_density.f08
    src/module_format_xsf_densities.f08
    src/module_check_noncollinear_XC_functional.f08
    src/module_gen_dens_grids_from_gaussian_basis_set_coeff.f08
    src/module_atomic_symbol_to_number.f08
    src/module_compute_dominant_atom_volumes.f08
    src/module_format_vasp_densities.f08
    src/module_run_valence_core_densities.f08
    src/module_core_iterator.f08
    src/module_local_multipole_moment_analysis.f08
    src/module_update_atomic_densities.f08
    src/module_DDEC3_valence_iterator.f08
    src/module_DDEC6_valence_iterator.f08
    src/module_compute_center_of_mass.f08
    src/module_total_multipole_moment_analysis.f08
    src/module_cloud_penetration.f08
    src/module_atomic_number_to_symbol.f08
    src/module_generate_atomic_polarizability_upper_bound_file.f08
    src/module_generate_net_atomic_charge_file.f08
    src/module_spin_functions.f08
    src/module_generate_atomic_spin_moment_file.f08
    src/module_collinear_spin_moments_iterator.f08
    src/module_noncollinear_spin_moments_iterator.f08
    src/module_generate_bond_order_file.f08
    src/module_determine_whether_pair_overlap_is_significant.f08
    src/module_calculate_final_BOs.f08
    src/module_integrate_bonding_terms.f08
    src/module_initialize_bond_pair_matrix.f08
    src/module_compute_local_atomic_exchange_vectors.f08
    src/module_prepare_BO_density_grids.f08
    src/module_perform_bond_order_analysis.f08
    src/module_compute_atomic_radial_moments.f08
    src/module_print_overlap_populations.f08
    src/module_print_atomic_densities_file.f08
    src/chargemol.f08
    )

set_source_files_properties(${SRC} PROPERTIES LANGUAGE Fortran)
set(CHARGEMOL_BINARY_NAME chargemol.x)

if(ENABLE_OPENMP)
    find_package(OpenMP REQUIRED)
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${OpenMP_Fortran_FLAGS}")
endif(ENABLE_OPENMP)

add_executable(${CHARGEMOL_BINARY_NAME} "${SRC}")

if(ENABLE_STATIC_LINK)
    # Only handle GCC for now.
    if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
        set(CHARGEMOL_LINK_FLAGS "-static-libgfortran")
    endif()
    if(DEFINED CHARGEMOL_LINK_FLAGS)
        set_target_properties(${CHARGEMOL_BINARY_NAME} PROPERTIES LINK_FLAGS ${CHARGEMOL_LINK_FLAGS})
    endif(DEFINED CHARGEMOL_LINK_FLAGS)
endif(ENABLE_STATIC_LINK)
