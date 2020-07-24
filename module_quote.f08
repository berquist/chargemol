 MODULE module_quote
 ! Copyright (c) 2014, 2015, 2016, 2017 by Thomas A. Manz and Nidia Gabaldon Limas. Rights reserved.
 
 USE module_common_variable_declarations

 IMPLICIT NONE
 
 CONTAINS
 SUBROUTINE quote()
 !===================================================================================

 WRITE(output_FID,*)'Normal termination of Chargemol program.'
 WRITE(output_FID,*)'Copyright (c) 2014, 2015, 2016, 2017 by Thomas A. Manz and Nidia Gabaldon Limas. Rights reserved.'
 WRITE(output_FID,*)'Use and distribution of this program is subject to certain licensing restrictions.'
 WRITE(output_FID,*)'Please see ddec.sourceforge.net for details.'
 WRITE(output_FID,*)' '
 WRITE(output_FID,*)'Please cite the following journal references when using this program to compute DDEC6 charges:'
 WRITE(output_FID,*)'N. Gabaldon Limas and T. A. Manz,  "Introducing DDEC6 atomic population analysis: part 2. '
 WRITE(output_FID,*)'Computed results for a wide range of periodic and nonperiodic materials," RSC Advances, Vol. 6 (2016) &
 &45727-45747.'
 WRITE(output_FID,*)'T. A. Manz and N. Gabaldon Limas, "Introducing DDEC6 atomic population analysis: part 1.'
 WRITE(output_FID,*)'Charge partitioning theory and methodology," RSC Advances, Vol. 6 (2016) 47771-47801.'
 WRITE(output_FID,*)'If desired, you can also cite references listed below for the DDEC3 and earlier methods.'
 WRITE(output_FID,*)' ' 
 WRITE(output_FID,*)'Please cite the following journal references when using this program to compute DDEC, Hirshfeld, or ISA charg&
 &es:'
 WRITE(output_FID,*)'Thomas A. Manz and David S. Sholl, "Improved Atoms-in-Molecule Charge Partitioning Functional for Simultaneou&
 &sly Reproducing the Electrostatic'
 WRITE(output_FID,*)'Potential and Chemical States in Periodic and Non-Periodic Materials", J. Chem. Theory Comput.&
 &, Vol. 8 (2012) 2844-2867.'
 WRITE(output_FID,*)'Thomas A. Manz and David S. Sholl, "Chemically Meaningful Atomic Charges that Reproduce the Electrostatic Pot&
 &ential in Periodic and Nonperiodic'
 WRITE(output_FID,*)'Materials", J. Chem. Theory Comput., Vol. 6 (2010) 2455-2468.'
 WRITE(output_FID,*)' '
 WRITE(output_FID,*)'Please cite the following journal reference when using this program to compute atomic spin moments:'
 WRITE(output_FID,*)'Thomas A. Manz and David S. Sholl, "Methods for Computing Accurate Atomic Spin Moments for Collinear and Nonc&
 &ollinear Magnetism in Periodic and'
 WRITE(output_FID,*)'Nonperiodic Materials", J. Chem. Theory Comput., Vol. 7 (2011) 4146-4164.'
 WRITE(output_FID,*)' '
 WRITE(output_FID,*)'Please cite the following reference when using this program to compute bond orders:'
 WRITE(output_FID,*)'Thomas A. Manz, "Introducing DDEC6 atomic population analysis: part 3. Comprehensive method to compute &
 &bond orders.", RSC Advances, Vol. 7 (2017) 45552-45581.'
 WRITE(output_FID,*)' '
 WRITE(output_FID,*)'Exiting Chargemol'

 END SUBROUTINE quote
 END MODULE module_quote