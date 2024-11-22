Getting Started
=====================
Sparcle_QC is a Python package that automatically creates QM/MM-like input files for the quantum chemistry codes Psi4, Q-Chem, and NWChem. It specifically is designed for use with protein:ligand complexes. Main features include cutting the QM region out of the protein, capping the QM region with hydrogen link atoms, obtaining point charges for the MM region, and altering QM/MM boundary charges.

Detailed below is the main workflow of Sparcle_QC with input keywords in parenthesis as they come up. These keywords must be provided at runtime to Sparcle_QC either by a dictionary in python or an input file on the command line as detailed in :doc:`installation`.

Sparcle_QC is fully automated. Users have the option of utilizing a CHARMM or Amber forcefield. As seen below, the Amber preparation can be done directly through Sparcle_QC, but CHARMM-GUI files must be uploaded manually to the working directly when using a CHARMM forcefield. 



1. System Prepartion
---------------------

Provided a protein-ligand complex (pdb_file), Sparcle_QC will create a new directory with the name of the input file to place all created files. 

If specified (pre-capped), when utilizing amber forcefields, Sparcle_QC can automatically cap the terminal residues with ACE and NME groups. This is not an option with CHARMM because users have an option to add capping groups when creating the input files in CHARMM-GUI. 

.. figure:: _static/light_termini.png
   :align: center
   :class: only-light
   :width: 600px
   
   In this example, the pre-capped keyword has been set to false so the termini are capped with ACE/NME. Here, the C terminus can be seen before and after capping with NME. 

.. figure:: _static/dark_termini.png
   :align: center
   :class: only-dark
   :width: 600px
   
   In this example, the pre-capped keyword has been set to false so the termini are capped with ACE/NME. Here, the C terminus can be seen before and after capping with NME. 


2. Cutting out the QM Region:
-----------------------------

To create a fully quantum mechanical subsystem, Sparcle_QC will carve out a chunk of the protein based on a user specified distance (cutoff_radius). This radius, in Angstroms, defines which atoms are in the QM region. If no specification is made, the radius will be measured from the closest atom in the ligand. However, if desired, a specific atom can be chosen by specifying the atomid (seed) and corresponding PDB (seed_file) to grow the QM region from. Sparcle_QC will make a cutout of the complex by severing C-C bonds, specifically the bond between the C of the peptide bond and the alpha C in the protein backbone.
   
.. figure:: _static/light_cutout.png
   :align: center
   :class: only-light
   :width: 600px
   
   The full system provided to Sparcle_QC can be seen on the left. After specifying a cutoff_radius of 3.5 Angstroms, the following cutout of the system was created. 

.. figure:: _static/dark_cutout.png
   :align: center
   :class: only-dark
   :width: 600px
   
   The full system provided to Sparcle_QC can be seen on the left. After specifying a cutoff_radius of 3.5 Angstroms, the following cutout of the system was created. 

Because C-C bonds have now been severed, some carbons on the exterior of the created cutout do not have four bonds and therefore have unsatisfied valencies. Sparcle_QC automatically caps these bonds with hydrogen link atoms. 

.. figure:: _static/light_cap.png
   :align: center
   :class: only-light
   :width: 600px
   
   In this zoomed in picture of the cutout from above, several carbonyls are picture with only three bonds because their fourth bonds were severed when creating the cutout. On the right, Sparcle_QC has capped these bonds with hydrogen link atoms. 

.. figure:: _static/dark_cap.png
   :align: center
   :class: only-dark
   :width: 600px
   
   In this zoomed in picture of the cutout from above, several carbonyls are picture with only three bonds because their fourth bonds were severed when creating the cutout. On the right, Sparcle_QC has capped these bonds with hydrogen link atoms. 

3. Obtaining Point Charges from a Forcefield
--------------------------------------------

The rest of the atoms, that were not included in the QM cutout, will be turned into external charges in the QM calculation so that they may partcipate in long range electrostatic interactions and polarize the orbitals in the QM region. 

.. figure:: _static/light_add_charges.png
   :class: only-light
   :align: center
   :width: 600px
   
   The same cutout from above, now capped with hydrogens, is picture on the left. After adding back the rest of the system as point charges, the system on the right is obtained.

.. figure:: _static/dark_add_charges.png
   :class: only-dark
   :align: center
   :width: 600px
   
   The same cutout from above, now capped with hydrogens, is picture on the left. After adding back the rest of the system as point charges, the system on the right is obtained.
   
Based on a user chosen forcefield, Sparcle_QC will automatically obtain the correct point charges for the system and add them to the QM cutout. As previously mentioned, Sparcle_QC is fully automated with amber forcefields (amber_ff) and will obtain the charges based on the specified forcefield, but if a CHARMM forcefield is requested, point charges can be obtained using CHARMM-GUI and uploaded with the same name as (pdb_file) but with the .psf extension. Along with the psf, users must provide paths to the CHARMM topology and parameter files (charmm_rtf and charmm_prm) and the ligand in its own file as ligand.pdb.  

For amber, a specific water model can be chosen (water_model) and/or the individual charges for the Oxygen, Hydrogen, and Extra Point Charge for a 4 point model can be chosen (o_charge, h_charge, ep_charge). 

4. Redistributing Boundary Charges
----------------------------------
   
Since the hydrogen link atoms were artifically added where there used to be a C-C bond, there are now point charges way too close to the QM region. This will cause overpolarization of the QM region. Sparcle_QC supports 9 different options to redistribute these boundary charges (charge_scheme). Read more about the different charge schemes here :doc:`user_guide`. 

.. figure:: _static/light_charges.png
   :align: center
   :class: only-light
   :width: 600px
   
   The capped subsystem along with the point charges is pictured on the left. This zoomed in view reveals that the artifical hydrogen link atom is several tenths of an angstrom from a point charge. On the right, Sparcle_QC has employed the Z1 charge scheme, zeroing the first MM atom next to the hydrogen, to avoid potential overpolarization. 

.. figure:: _static/dark_charges.png
   :align: center
   :class: only-dark
   :width: 600px
   
   The capped subsystem along with the point charges is pictured on the left. This zoomed in view reveals that the artifical hydrogen link atom is several tenths of an angstrom from a point charge. On the right, Sparcle_QC has employed the Z1 charge scheme, zeroing the first MM atom next to the hydrogen, to avoid potential overpolarization. 


5. Writing a QM Input File
---------------------------

Sparcle_QC will write input files for 3 different softwares: NWChem, Psi4, Q-CHEM (software) for a user chosen level of theory (method and basis_set). To avoid ambiguities in ligand isomers, the charge of the ligand must be specified (ligand_charge). Users may also specify if ghost atoms should be added for counterpoise correction (cp). 

Software specific keywords may also be controlled via the Sparcle_QC input file. More details can be found here: :doc:`api`.

6. Extra Features
-----------------

Extra amber forcefields (other_amber_ff) needed to obtain proper point charges for the protein can be added as well, for example for phosphorylated amino acids. 

If using FI-SAPT in Psi4, Sparcle_QC can automatically create the functional group partition files (fisapt_partition), fA.dat and fB.dat. More details about FI-SAPT can be found on `Psi4's website <https://psicode.org/psi4manual/master/fisapt.html>`_

.. note::
   This partition has been successful in our group for the systems we have come across, but it is important to check these partition files for system specific nuances of functional groups. 

.. figure:: _static/fsapt_sparcle.png
   :align: center
   :class: only-light
   :width: 600px
   
   Users can utilize FI-SAPT in Psi4 to obtain functional group partitions of the interaction energy. After running Psi4's post processing script for fsapt, a heat map of the interaction energy can be visualized on the protein structure, where red indicates a favorable interaction and blue indicates an unfavorable interaction.

.. figure:: _static/dark_fsapt_sparcle.png
   :align: center
   :class: only-dark
   :width: 600px
   
   Users can utilize FI-SAPT in Psi4 to obtain functional group partitions of the interaction energy. After running Psi4's post processing script for fsapt, a heat map of the interaction energy can be visualized on the protein structure, where red indicates a favorable interaction and blue indicates an unfavorable interaction.

Sparcle_QC can also match QM regions between two seperate runs. Users can provide a previous Sparcle_QC runs finalized cx_autocap_fixed.pdb (template_path) to serve as a template for the current Sparcle_QC run. The residues of the current protein will be mapped to the residues of the template protein by their neighborhood and the same amount of each residue will be included in the QM regions. Anything that cannot be matched via the residue neighborhood between the two proteins will be treated as normal with a distance based criteria. Details can be found here :doc:`api`.

.. note::
   This conversion script has been successful for system in our group. However, waters and uncommon atoms between the systems will purely be treated with the distance cutoff, so it is important to ensure that the CAPPED_qm.pdb files match in the way that the user intended after Sparcle_QC completes. 

7. Syntax of the Input Files
-----------------------------
These keywords that have been specified here in parenthesis should be assembled in an input file with any extension with the keyword followed by a colon and the value. Example: 
        
        amberff: ff19SB

Full example input files can be found here: :doc:`api`

