User Guide
===============

General Information
*********************

.. dropdown:: Cutting and Capping Bonds
    
    Alpha carbon - carbonyl carbon bonds are cut to separate the QM region from the MM region. 
   
    .. image:: _static/cut-lightmode.svg
        :class: only-light
        :align: center
        :width: 500 px

    .. image:: _static/cut-darkmode.svg
        :class: only-dark
        :align: center
        :width: 500 px

    The valencies on the QM atoms of the cut bonds are satisfied by placing a hydrogen along the cut bond. The hydrogen bond length, :math:`r_{\rm{Q1-H1}}`, is determined by

    .. math::
    
        r_{\rm{Q1-H1}} = r_{\rm{Q1-M1}}\frac{r^{\rm{QM0}}_{\rm{Q1-HL}}}{r^{\rm{MM0}}_{\rm{Q1-M1}}},
    
    where :math:`r_{\rm{Q1-M1}}` is the bond length of the cut bond, :math:`r^{\rm{QM0}}_{\rm{Q1-M1}}` is the corresponding bond length according to the force field used, and :math:`r^{\rm{QM0}}_{\rm{Q1-HL}}` is the force field bond length for the hydrogen link bond. [1]_

    .. [1]      L.H. Hu, P. Soderhjelm, and U. Ryde, "On the convergence of QM/MM energies," *J. Chem. Theory Comput.* 7, 761-777 (2011). `<https://doi.org/10.1021/ct100530r>`_

.. dropdown:: QM/MM Boundary Charges

    Point charges too close to the capped QM region may cause overpolarization. Users have the option of choosing one of nine charge schemes to alter charges at this boundary. These schemes are:

    * **Z1**: Charges of the M1 atoms are set to zero.
    * **Z2**: Charges of the M1 and M2 atoms are set to zero.
    * **Z3**: Charges of the M1, M2, and M3 atoms are set to zero.
    * **DZ1**: At each cut, the charge of M1 is set to zero, and the charge needed to return the MM boundary residue to its original integer charge is evenly distributed to all MM atoms in that residue.
    * **DZ2**: At each cut, the charges of M1 and M2 atoms are set to zero, and the charge needed to return the MM boundary residue to its original integer charge is evenly distributed to all MM atoms in that residue.
    * **DZ3**: At each cut, the charges of M1, M2, and M3 atoms are set to zero, and the charge needed to return the MM boundary residue to its original integer charge is evenly distributed to all MM atoms in that residue.
    * **BRC**: At each cut, the charge of M1 is set to zero, and the charge needed to return the MM boundary residue to its original integer charge is evenly distributed to the midpoints of the M1-M2 bonds.
    * **BRCD**: At each cut, the charge of M1 is set to zero, and the charge needed to return the MM boundary residue to its original integer charge is evenly distributed to the midpoints of the M1-M2 bonds, but doubled. This charge is also subtracted from each M2 atom within the residue.
    * **BRC2**: At each cut, the charge of M1 is set to zero, and the charge needed to return the MM boundary residue to its original integer charge is evenly distributed to the M2 atoms in that residue.

    For SAPT0 in Psi4, we recommend BRC. [2]_ 

    .. image:: _static/qmmm-lightmode.svg
        :class: only-light
        :align: center
        :width: 700 px
        
    .. image:: _static/qmmm-darkmode.svg
        :class: only-dark
        :align: center
        :width: 700 px

    .. [2] C.S. Glick, A. Alenaizan, D.L. Cheney, C.E. Cavender, and C.D. Sherrill, "Electrostatically embedded symmetry-adapted perturbation theory," *J. Chem. Phys.* 161, 134112 (2024). `<https://doi.org/10.1063/5.0221974>`_

.. dropdown:: Resulting Input File for QM Calculation

    The result of running Sparcle_QC is an input file for either Psi4, Q-Chem, or NWChem. In general, the input file will look similar to the files below.

    .. tab-set::
    
        .. tab-item:: Psi4
    
            .. code-block:: 
                
                """
                This Psi4 file was created using Sparcle-QC with the following specifications:
                pdb_file: 2cji.pdb
                cutoff: 5
                ... 
                ... # copy of Sparcle_QC input file
                """
                
                import psi4
                import numpy as np
                import qcelemental as qcel
                import time
                
                start = time.time()
                
                psi4.set_memory('60 GB')
                psi4.core.set_num_threads(2)
                
                psi4.core.set_output_file('psi4_file.out', False)
                dimer = psi4.geometry('''
                0 1
                 N -17.183 -79.238 -85.266
                 C -13.352 -80.694 -86.001
                 O -17.152 -80.421 -85.511
                 H -14.073 -80.711 -88.421
                 H -12.998 -81.949 -89.060
                 H -8.563 -81.173 -79.793
                 H -7.409 -80.135 -79.552
                --
                1 1
                 C -17.273 -84.206 -80.622
                 O -16.663 -84.413 -79.570
                 C -16.682 -81.881 -81.407
                 C -16.314 -82.218 -82.856
                 C -17.017 -80.384 -81.352
                 H -18.036 -81.829 -79.079
                 H -18.516 -82.920 -81.796
                 H -15.811 -82.046 -80.774
                 H -15.543 -81.535 -83.213
                 H -15.939 -83.235 -82.932
                 H -17.198 -82.118 -83.486
                units angstrom
                symmetry c1
                no_com
                no_reorient
                ''')
                
                Chargefield_B = np.array([
                0.5972,-25.097,-92.541,-80.98
                ,-0.5679,-26.081,-91.792,-81.032
                ,-0.3662,-24.383,-92.801,-79.671
                ,0.1123,-25.065,-93.243,-78.959
                ,0.1123,-23.555,-93.476,-79.829
                ,0.1123,-24.006,-91.874,-79.266]).reshape((-1,4))
                Chargefield_B[:,[1,2,3]] /= qcel.constants.bohr2angstroms
                
                psi4.set_options({
                'basis': 'aug-cc-pv(D+d)z',
                'freeze_core':'true',
                'scf_type':'df'
                })
                
                e = psi4.energy('sapt0', external_potentials={'B':Chargefield_B})
                
                end=time.time()
                wall_time = '{:.2f}'.format(float(end-start))
                with open ('psi4_file.out', 'a') as output:
                    output.write(f'Wall time: {wall_time} seconds')

    
        .. tab-item:: Q-Chem
    
            .. code-block:: 

                """
                This Q-Chem file was created using Sparcle-QC with the following specifications:
                pdb_file: 2cji.pdb
                cutoff: 5
                ...
                ... # copy of Sparcle_QC input file
                """
                
                $molcule
                4 1
                 N -17.183 -79.238 -85.266
                 C -13.352 -80.694 -86.001
                 O -17.152 -80.421 -85.511
                 H -14.073 -80.711 -88.421
                 H -12.998 -81.949 -89.060
                 H -8.563 -81.173 -79.793
                 H -7.409 -80.135 -79.552
                 C -17.408 -77.515 -77.251
                 O -16.597 -76.592 -77.177
                 N -18.231 -77.603 -78.308
                 C -18.398 -76.535 -79.306
                 C -19.485 -75.554 -78.882
                 O -20.421 -75.955 -78.193
                 H -18.875 -78.381 -78.368
                 H -17.467 -75.980 -79.429
                 H -18.673 -76.965 -80.270
                 N -19.491 -74.326 -79.419
                $end
                
                $external_charges
                    -25.097    -92.541    -80.98    0.5972
                    -26.081    -91.792    -81.032    -0.5679
                    -24.383    -92.801    -79.671    -0.3662
                    -25.065    -93.243    -78.959    0.1123
                    -23.555    -93.476    -79.829    0.1123
                    -24.006    -91.874    -79.266    0.1123
                $end
                
                $rem
                METHOD hf
                BASIS 6-31g*
                JOBTYPE sp
                $end
    
    
        .. tab-item:: NWChem
    
            .. code-block:: 
    
                """
                This NWChem file was created using Sparcle-QC with the following specifications:
                pdb_file: 2cji.pdb
                ...
                ... # copy of Sparcle_QC input file
                """
                
                START
                SCRATCH_DIR /scratch/user/
                PERMANENT_DIR /scratch/user/
                MEMORY 32 GB
                
                geometry nocenter noautoz noautosym
                4 1
                 N -17.183 -79.238 -85.266
                 C -13.352 -80.694 -86.001
                 O -17.152 -80.421 -85.511
                 H -14.073 -80.711 -88.421
                 H -12.998 -81.949 -89.060
                 H -8.563 -81.173 -79.793
                 H -7.409 -80.135 -79.552
                 C -17.408 -77.515 -77.251
                 O -16.597 -76.592 -77.177
                 N -18.231 -77.603 -78.308
                 C -18.398 -76.535 -79.306
                 C -19.485 -75.554 -78.882
                 O -20.421 -75.955 -78.193
                 H -18.875 -78.381 -78.368
                 H -17.467 -75.980 -79.429
                 H -18.673 -76.965 -80.270
                 N -19.491 -74.326 -79.419
                end
                
                bq
                    -25.097    -92.541    -80.98    0.5972
                    -26.081    -91.792    -81.032    -0.5679
                    -24.383    -92.801    -79.671    -0.3662
                    -25.065    -93.243    -78.959    0.1123
                    -23.555    -93.476    -79.829    0.1123
                    -24.006    -91.874    -79.266    0.1123
                end
                
                basis
                * library cc-pvdz
                end
                
                task hf energy


Options
********

.. dropdown:: Required
    
    .. csv-table:: 
        :file: options/required.csv

.. dropdown:: Force Fields

    .. csv-table:: 
        :file: options/forcefields.csv

.. dropdown:: Psi4

    .. csv-table:: 
        :file: options/psi4.csv

.. dropdown:: Q-Chem

    .. csv-table:: 
        :file: options/qchem.csv

.. dropdown:: NWChem

    .. csv-table:: 
        :file: options/nwchem.csv

.. dropdown:: Other

    .. csv-table:: 
        :file: options/other.csv


Example Inputs
***************

.. dropdown:: F-SAPT with Psi4 and Amber

    The following input file will create an F-SAPT file for the protein:ligand complex to be run with Psi4. It will also create the functional group partitions needed for post-processing, fA.dat and fB.dat. 
    
    .. code-block::
 
        pdb_file: 2cji.pdb
        pre-capped: true
        cutoff: 5
        seed: ligand
        charge_scheme: BRC
        ligand_charge: 0
        method: fisapt0
        fisapt_partition: true
        basis_set: aug-cc-pv(D+d)z
        amber_ff: ff19SB
        water_model: opc
        o_charge: 0
        h_charge: 0.6791
        ep_charge: -1.3582
        software: psi4
        mem: 60 GB
        nthreads: 10


.. dropdown:: B3LYP with Q-Chem and CHARMM

    The following input file will create prepare 3 Q-Chem files: one with the ligand  (fully QM), one with the protein (QM/MM), and one with the complex (QM/MM). These could be used to calculate a supermolecular interaction energy. We will turn on counterpoise correction, which will include ghost atoms for the QM dimer in all three files.
    
    .. code-block::
 
        pdb_file: 3qxp.pdb
        cutoff: 5
        seed: ligand
        charge_scheme: DZ3
        ligand_charge: 0
        method: b3lyp
        basis_set: 6-31G*
        charmm_rtf: top_all36_prot.rtf
        charmm_prm: par_all36m_prot.prm
        water_model: tip3p
        software: q-chem

.. dropdown:: HF with NWChem and Amber

    The following input file will create prepare 3 NWChem files: one with the ligand  (fully QM), one with the protein (QM/MM), and one with the complex (QM/MM). These could be used to calculate a supermolecular interaction energy. We will turn on counterpoise correction, which will include ghost atoms for the QM dimer in all three files. The QM region will grow starting from a single ligand atom.
    
    .. code-block::

        pdb_file: 2cji.pdb
        pre-capped: true
        cutoff: 8.5
        seed: 4247
        seed_file: 4yff.pdb
        charge_scheme: BRC
        ligand_charge: 0
        method: hf 
        basis_set: aug-cc-pv(D+d)z
        amber_ff: ff19SB
        water_model: opc
        o_charge: 0
        h_charge: 0.6791
        ep_charge: -1.3582
        software: nwchem
        nwchem_scratch: /scratch/user
        nwchem_perm: /scratch/user
        mem: 60 GB

.. dropdown:: Templating a QM Region for Congeneric Ligands

    Studies that compare a protein with two similar ligand structures may choose to equilibrate protein structures for each ligand. In this case, the two PDBs may be similar in structure, but not identical, and their coordinates likely will not match. Here, we show the steps of (1) creating a SAPT input file for one ligand (named methyl), then (2) using the QM region of 'methyl' as a template for cutting the QM region of the other ligand, named 'chlorine'.


    Step 1, the following is 'methyl.in':

    .. code-block::

        pdb_file: 2cji_methyl.pdb
        pre-capped: true
        cutoff: 5
        seed: ligand
        charge_scheme: BRC
        ligand_charge: 0
        method: fisapt0
        basis_set: aug-cc-pv(D+d)z
        amber_ff: ff19SB
        water_model: opc
        o_charge: 0
        h_charge: 0.6791
        ep_charge: -1.3582
        software: psi4
        mem: 60 GB
        nthreads: 10
 
       
    Step 2, the following is 'chlorine.in':

    .. code-block::

        pdb_file: 2cji_chlorine.pdb
        pre-capped: true
        template_path: methyl/cx_autocap_fixed.pdb
        charge_scheme: BRC
        ligand_charge: 0
        method: fisapt0
        basis_set: aug-cc-pv(D+d)z
        amber_ff: ff19SB
        water_model: opc
        o_charge: 0
        h_charge: 0.6791
        ep_charge: -1.3582
        software: psi4
        mem: 60 GB
        nthreads: 10
  
    With the two SAPT files, a relative interaction energy can be computed, giving insight into which ligand is more stable within the protein pocket. 
      
.. dropdown:: Convergence Study with Increasing QM Region Size via API

    A Python loop can be used to generate multiple input files with an increasing size of the QM region. We can increase the size of the QM region by incrementing the cutoff. An example Python script is below.    

    .. code-block::
 
        import sparcle_qc

	inputs = {
            'pdb_file': '2cji.pdb',
            'pre-capped': 'True',
            'seed': 'ligand',
            'charge_scheme': 'BRC',
            'ligand_charge': 0,
            'method': 'fisapt0',
            'basis_set': 'aug-cc-pv(D+d)z',
            'amber_ff': 'ff19SB',
            'water_model': 'opc' ,
            'o_charge': 0,
            'h_charge': 0.6791,
            'ep_charge': -1.3582,
            'software': 'psi4',
            'mem': '60 GB',
            'nthreads': 10}

	cutoffs = [3, 4, 5]

	for c in cutoffs:
	    inputs['cutoff'] = f'{c}'
	    inputs['input_filename'] = f'cutoff_{c}.in'
	    print(inputs)
	    sparcle_qc.run_sparcle(user_options = inputs)

                                                     
