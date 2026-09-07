This README.txt file was generated on <20260720> by Ilkwon Cho.

#
# General instructions for completing README:
# Please leave all commented sections in README (do not delete any text).
#

-------------------
GENERAL INFORMATION
-------------------

1. Title of Dataset:

Benchmark datasets for label-free, quantum-mechanically informed scoring of
protein-ligand complexes with machine-learned interatomic potentials

#
# Authors: Include contact information for at least the first author and
# corresponding author (if not the same), specifically email address, phone
# number (optional, but preferred), and institution. Contact information for
# all authors is preferred.
#

2. Author Information

First Author Contact Information
    Name: Ilkwon Cho
    ORCiD: 0000-0002-7279-7384
    Institution: Carnegie Mellon University, Department of Chemistry
    Address: 5000 Forbes Avenue, Pittsburgh, PA 15213, USA
    Email: ilkwonc@andrew.cmu.edu

Author Contact Information
    Name: Hatice Gokcan
    ORCiD: 0000-0002-0112-5679
    Institution: Carnegie Mellon University, Department of Chemistry
    Address: 5000 Forbes Avenue, Pittsburgh, PA 15213, USA
    Email: hgokcan@andrew.cmu.edu 

Corresponding Author Contact Information
    Name: Olexandr Isayev
    ORCiD: 0000-0001-7581-8497
    Institution: Carnegie Mellon University, Departments of Chemistry and
                 Materials Science and Engineering
    Address: 5000 Forbes Avenue, Pittsburgh, PA 15213, USA
    Email: olexandr@olexandrisayev.com

---------------------
DATA & FILE OVERVIEW
---------------------

#
# Directory of Files in Dataset: List and define the different files included
# in the dataset. This serves as its table of contents.
#

Directory of Files:

   A. Directory: KIN66/
      Short description: A 66-system kinase protein-ligand interaction-energy
      benchmark containing CDK2 (16 systems), P38 (34 systems), and TYK2
      (16 systems). The directory contains prepared complex, protein-fragment,
      and ligand geometries, together with CSV tables of DFT, machine-learned
      interatomic potential, and semiempirical interaction energies.

      A.1 Filename: KIN66/cdk2.csv
          Short description: Interaction energies for the 16 CDK2 systems.

      A.2 Filename: KIN66/p38.csv
          Short description: Interaction energies for the 34 P38 systems.

      A.3 Filename: KIN66/tyk2.csv
          Short description: Interaction energies for the 16 TYK2 systems.

      A.4 Directories: KIN66/cdk2/, KIN66/p38/, KIN66/tyk2/
          Short description: Per-system structures in PDB, SDF, and XYZ
          formats. For every molecule identifier, the complex, isolated
          protein fragment, and isolated ligand preserve the coordinates used
          to calculate the supramolecular interaction energy.

   B. Directory: FEP+_R-group_set/
      Short description: Scoring results and prepared structures for eight
      congeneric FEP+ R-group targets: BACE, CDK2, JNK1, MCL1, P38, PTP1B,
      thrombin, and TYK2. Experimentally assayed benchmark ligands are labeled
      as actives. The compounds labeled inactive are target-specific putative
      decoys generated with DeepCoy; they are not experimentally confirmed
      nonbinders.

      B.1 Filenames: <target>_active.csv
          Short description: Results for the experimentally annotated subset
          used for affinity-correlation analysis. Across all eight files this
          subset contains 199 ligands. The target counts are BACE 36, CDK2 16,
          JNK1 21, MCL1 42, P38 34, PTP1B 23, thrombin 11, and TYK2 16.

      B.2 Filenames: <target>_all.csv
          Short description: Combined active and putative-decoy tables used
          for pose-conditioned virtual-screening/enrichment analysis. These
          files contain 5,970 DeepCoy decoys, corresponding to 30 decoys for
          each of the 199 annotated actives.

      B.3 Filename: geometry.tar.gz
          Short description: Gzip-compressed tar archive containing the
          prepared active and putative-decoy structures. It can be extracted
          with:

              tar -xzf geometry.tar.gz

          Within the archive, structures are organized as
          geometry/{actives,inactives}/<target>/. Each system is represented
          by a ligand SDF and a corresponding truncated protein PDB. The
          ligand and protein coordinates are in the same coordinate frame and
          together define the scored protein-ligand complex.

   C. Directory: PLA15/
      Short description: Results for the external 15-system PLA15 charged
      protein-ligand benchmark. The original PLA15 geometries were not created
      by the authors of the present dataset. New full-system single-point DFT
      interaction energies and model predictions were evaluated at the
      published geometries.

      C.1 Filename: PLA15/PLA15.csv
          Short description: Original PLA15 reference interaction energies,
          new DFT interaction energies, and model predictions for 15 systems.

      C.2 Directory: PLA15/geometry/
          Short description: PLA15 complex and fragment geometries in PDB and
          XYZ formats. Original PLA15 resource and provenance:

              http://cuby4.molecular.cz/dataset_pla15.html

          Users of these geometries should acknowledge and cite the original
          PLA15 resource and publication.

   D. Directory: HiQBind/
      Short description: A curated 108-complex subset of the HiQBind
      noncovalent protein-ligand dataset used to evaluate affinity ranking on
      experimentally resolved co-crystal poses. The subset contains seven
      targets: trmD (25 complexes), HRAS (13), HSP90AA1 (22), dsbA (12),
      RXRA (10), WDR5 (13), and CREBBP (13).

      D.1 Filename: HiQBind/hiqbind.csv
          Short description: Experimental dissociation constants, ligand and
          protein identifiers, activity labels, AIMNet2 metrics, and baseline
          scoring results for the 108 co-crystal complexes.

      D.2 Filename: HiQBind/geometry.tar.gz
          Short description: Gzip-compressed tar archive containing the
          target-organized co-crystal ligand poses and corresponding
          truncated, hydrogen-capped protein fragments. It can be extracted
          with:

              tar -xzf geometry.tar.gz

          Within the archive, structures are organized as
          geometry/{trmD,HRAS,HSP90AA1,dsbA,RXRA,WDR5,CREBBP}/. Each PDBID is
          represented by one ligand SDF and one protein PDB in the same
          coordinate frame. The 108 CSV rows, ligand files, and protein files
          are in exact one-to-one correspondence.

          This subset was derived from the HiQBind resource described in:

              Wang, Y. et al. A workflow to create a high-quality
              protein-ligand binding dataset for training, validation, and
              prediction tasks. Digital Discovery 4, 1209-1220 (2025).
              https://doi.org/10.1039/D4DD00357H

          Users should cite the original HiQBind resource when reusing these
          structures or associated experimental binding data.


# Additional Notes on File Relationships, Context, or Content
# (for example, if a user wants to reuse and/or cite your data, what
# information would you want them to know?):

Interaction energies were evaluated with the supramolecular expression

    E_interaction = E_complex - E_protein - E_ligand

using complex and isolated-fragment single-point energies at identical atomic
coordinates. Unless otherwise noted below, interaction-energy and binding-
energy columns are reported in kcal mol^-1.

KIN66 and PLA15 provide DFT interaction-energy benchmarks. The FEP+ R-group
set provides affinity/scoring tables and prepared structures for correlation
and retrospective, pose-conditioned enrichment analyses. HiQBind provides an
independent affinity/scoring benchmark based on experimentally resolved
co-crystal poses. The FEP+ putative decoys should not be interpreted as
experimentally validated inactive compounds.

Blank CSV fields are missing values (NA). In FEP+_R-group_set, E_exp and FEP
are unavailable for generated decoys. Some AIMNet2 values are also blank when
structure validation or scoring did not produce a usable result. In KIN66,
directly converged omegaB97M-V values and linearly estimated values are stored
in separate columns.


#
# File Naming Convention: Define your File Naming Convention (FNC), the
# framework used for naming your files systematically to describe what they
# contain, which could be combined with the Directory of Files.
#

File Naming Convention:

KIN66

    <target>.csv
        Target-level interaction-energy table.

    <target>_complex_<molecule_id>.pdb
    <target>_complex_<molecule_id>.xyz
        Prepared truncated protein-ligand complex.

    <target>_protein_<molecule_id>.pdb
    <target>_protein_<molecule_id>.xyz
        Isolated capped protein fragment in its complex geometry.

    <target>_ligands_<molecule_id>.sdf
    <target>_ligands_<molecule_id>.xyz
        Isolated ligand in its complex geometry. The plural word "ligands"
        is retained in the filename even though each file contains one
        ligand.

FEP+_R-group_set

    <target>_active.csv
        Experimentally annotated affinity-correlation subset.

    <target>_all.csv
        Combined active and putative-decoy virtual-screening table.

    geometry/actives/<target>/<target>_ligands_<ligand_id>.sdf
    geometry/actives/<target>/<target>_protein_<ligand_id>.pdb
        Prepared active ligand and corresponding truncated protein fragment.

    geometry/inactives/<target>/<target>_ligands_decoy<number>.sdf
    geometry/inactives/<target>/<target>_protein_decoy<number>.pdb
        Prepared putative-decoy ligand and corresponding truncated protein
        fragment.

HiQBind

    hiqbind.csv
        Experimental affinities, system metadata, and scoring results for the
        curated 108-complex subset.

    geometry.tar.gz
        Compressed archive of all HiQBind subset geometries. Extract with:

            tar -xzf geometry.tar.gz

    Files within the extracted archive:

    geometry/<target>/<PDBID>_<ligand_name>_ligands.sdf
        Co-crystallized ligand coordinates.

    geometry/<target>/<PDBID>_<ligand_name>_protein.pdb
        Corresponding truncated and hydrogen-capped protein fragment in the
        same coordinate frame as the ligand SDF.

    Target directory names are display names mapped from UniProt identifiers:

        B1MDI3,B1MDI3  -> trmD
        P01112,Q07889   -> HRAS
        P07900          -> HSP90AA1
        P0AEG4,P0AEG4  -> dsbA
        P19793          -> RXRA
        P61964          -> WDR5
        Q92793          -> CREBBP

PLA15

    <system_id>.pdb
        Combined PLA15 complex. REMARK records specify the total charge,
        fragment charges, and selections used to separate fragments.

    <system_id>.xyz
        Combined complex coordinates.

    <system_id>_a.xyz
        Protein-fragment coordinates.

    <system_id>_b.xyz
        Ligand coordinates.


#
# Data Description: A data description, dictionary, or codebook defines the
# variables and abbreviations used in a dataset. This information can be
# included in the README file, in a separate file, or as part of the data
# file. If it is in a separate file or in the data file, explain where this
# information is located and ensure that it is accessible without specialized
# software. (We recommend using plain text files or tabular plain text CSV
# files exported from spreadsheet software.)
#

---------------------------------------
DATA DESCRIPTION FOR: KIN66/*.csv
---------------------------------------

1. Number of variables:

   cdk2.csv: 13
   p38.csv: 14
   tyk2.csv: 14

2. Number of cases/rows:

   cdk2.csv: 16
   p38.csv: 34
   tyk2.csv: 16
   Total: 66

3. Missing data codes:

   Blank field / NA: value not available. Directly converged
   int_wb97m_dzvpd values are available for all 16 CDK2 systems, one P38
   system, and three TYK2 systems. Estimated values for P38 and TYK2 are
   provided separately in int_wb97m_dzvpd_estimated.

4. Variable List

    A. Name: molecule_id
       Description: Target-specific molecule identifier used to match the CSV
       row to its geometry files.

    B. Names: MACE_OFF23_interaction, MACE_OFF24_interaction,
              MACE_OMol_interaction, MACE_Polar_L_interaction,
              UMA_S_interaction, esen_sm_conserving_interaction,
              gxtb_interaction, pm6_ml_interaction,
              AIMNet2_2025_interaction
       Description: Predicted protein-ligand interaction energies from the
       named MLIP or semiempirical model, in kcal mol^-1.

    C. Name: int_b973c
       Description: B97-3c DFT interaction energy, in kcal mol^-1.

    D. Name: int_wb97m
       Description: omegaB97M-D3(BJ)/def2-TZVPP DFT interaction energy, in
       kcal mol^-1.

    E. Name: int_wb97m_dzvpd
       Description: Directly converged omegaB97M-V/def2-TZVPD DFT interaction
       energy, in kcal mol^-1.

    F. Name: int_wb97m_dzvpd_estimated
       Description: Estimated omegaB97M-V/def2-TZVPD interaction energy for
       P38 and TYK2, obtained from the linear relationship with
       omegaB97M-D3(BJ)/def2-TZVPP. This column is not present in cdk2.csv,
       where all omegaB97M-V calculations converged directly.


----------------------------------------------------
DATA DESCRIPTION FOR: FEP+_R-group_set/*.csv
----------------------------------------------------

1. Number of variables:

   <target>_active.csv: 10
   <target>_all.csv: 11

2. Number of cases/rows:

   Target       active.csv      all.csv
   BACE                 36         1116
   CDK2                 16          496
   JNK1                 21          651
   MCL1                 42         1302
   P38                  34         1054
   PTP1B                23          713
   Thrombin             11          341
   TYK2                 16          496
   Total               199         6169

3. Missing data codes:

   Blank field / NA: not applicable or not available. E_exp and FEP are not
   assigned to putative decoys. Missing AIMNet2 values indicate that no usable
   score was produced, including cases rejected by post-minimization structural
   validation.

4. Variable List

    A. Name: Ligand
       Description: Target-specific ligand identifier or generated decoy
       identifier.

    B. Names: gnina, vina, smina, unidock
       Description: Raw scores produced by the named docking/scoring methods.
       GNINA is stored as its affinity output; Vina, Smina, and Uni-Dock are
       stored in their native score convention.

    C. Name: rtmscore
       Description: Raw RTMScore output. Present only in <target>_all.csv.

    D. Name: active
       Description: Boolean class label. True denotes a benchmark ligand;
       False denotes a DeepCoy-generated putative decoy. False does not mean
       experimentally confirmed inactivity.

    E. Name: AIMNet2(Eint)
       Description: Gas-phase protein-ligand interaction energy from
       AIMNet2(2025), in kcal mol^-1. Lower values indicate more favorable
       interactions.

    F. Name: AIMNet2(Score)
       Description: Composite endpoint score equal to interaction energy plus
       ligand desolvation and local conformational-strain contributions, in
       kcal mol^-1. Lower values indicate more favorable scores.

    G. Name: FEP
       Description: Reported FEP+ binding free-energy estimate, in
       kcal mol^-1.

    H. Name: E_exp
       Description: Experimental binding free energy, in kcal mol^-1.


------------------------------------------
DATA DESCRIPTION FOR: HiQBind/hiqbind.csv
------------------------------------------

1. Number of variables: 21

2. Number of cases/rows: 108

   Target       Number of complexes
   trmD                          25
   HRAS                          13
   HSP90AA1                      22
   dsbA                          12
   RXRA                          10
   WDR5                          13
   CREBBP                        13
   Total                        108

3. Missing data codes:

   Blank field / NA: value not available. boltz_kcalmol is unavailable for
   one complex. All other fields are complete in this release.

4. Variable List

    A. Name: PDBID
       Description: Four-character Protein Data Bank identifier used to match
       each row to its ligand and protein geometry files.

    B. Names: Ligand Name, Ligand Chain, Ligand Residue Number
       Description: PDB chemical-component identifier, chain identifier, and
       residue number of the co-crystallized ligand.

    C. Name: Binding Affinity Measurement
       Description: Experimental measurement type. All 108 entries are
       dissociation constants (kd).

    D. Names: Protein UniProtID, Protein UniProtName
       Description: UniProt accession identifier(s) and protein name supplied
       with the curated HiQBind record.

    E. Name: target
       Description: Short display name used in the accompanying manuscript
       and geometry directory structure: trmD, HRAS, HSP90AA1, dsbA, RXRA,
       WDR5, or CREBBP.

    F. Name: Ligand SMILES
       Description: SMILES representation of the ligand.

    G. Name: kd_nM
       Description: Experimental dissociation constant in nM.

    H. Name: Activity
       Description: Binary label derived from kd_nM. "active" denotes
       kd <= 10,000 nM; "Inactive" denotes kd > 10,000 nM. This subset
       contains 68 active and 40 inactive entries.

    I. Name: Log Binding Affinity nM
       Description: Natural logarithm of kd_nM, ln(kd_nM).

    J. Name: AIMNet2(Eint)
       Description: Gas-phase protein-ligand interaction energy from
       AIMNet2(2025), in kcal mol^-1. Lower values indicate more favorable
       interactions.

    K. Name: AIMNet2(Score)
       Description: Composite endpoint score equal to interaction energy plus
       ligand desolvation and local conformational-strain contributions, in
       kcal mol^-1. Lower values indicate more favorable scores.

    L. Name: Kd_kcalmol
       Description: Experimental dissociation constant converted to a binding
       free-energy convention in kcal mol^-1.

    M. Names: gnina, vina, smina, unidock
       Description: Raw scores produced by the named docking/scoring methods
       when rescoring the curated co-crystal geometries. These methods were
       not used to generate the HiQBind poses.

    N. Name: rtmscore
       Description: Raw RTMScore output for the co-crystal complex.

    O. Name: boltz_kcalmol
       Description: Boltz-2 affinity prediction converted to kcal mol^-1.


---------------------------------------
DATA DESCRIPTION FOR: PLA15/PLA15.csv
---------------------------------------

1. Number of variables: 12

2. Number of cases/rows: 15

3. Missing data codes: N/A

4. Variable List

    A. Name: molecule_id
       Description: PLA15 system identifier matching the geometry filenames.

    B. Name: PLA15_reference_interaction
       Description: Original PLA15 composite benchmark interaction energy, in
       kcal mol^-1.

    C. Names: MACE_OMol_interaction, MACE_Polar_L_interaction,
              UMA_S_interaction, esen_sm_conserving_interaction,
              gxtb_interaction, pm6_ml_interaction,
              AIMNet2_2025_interaction
       Description: Predicted protein-ligand interaction energies from the
       named MLIP or semiempirical model, in kcal mol^-1.

    D. Name: int_b973c
       Description: Full-system B97-3c DFT interaction energy evaluated at the
       original PLA15 geometry, in kcal mol^-1.

    E. Name: int_wb97m
       Description: Full-system omegaB97M-D3(BJ)/def2-TZVPP DFT interaction
       energy evaluated at the original PLA15 geometry, in kcal mol^-1.

    F. Name: int_wb97m_dzvpd
       Description: Full-system omegaB97M-V/def2-TZVPD DFT interaction energy
       evaluated at the original PLA15 geometry, in kcal mol^-1.


--------------------------
METHODOLOGICAL INFORMATION
--------------------------

#
# Software: If specialized software(s) generated your data or are necessary
# to interpret it, provide for each (if applicable): software name, version,
# system requirements, developer, and URLs. If you developed the software,
# provide the relevant source components or repository pointer.
#

1. Software-specific information:

Name: ORCA
Version: 6
System Requirements: N/A
Open Source? (Y/N): N
Executable URL: https://www.faccts.de/orca
Source Repository URL: N/A
Developer: FACCTs GmbH
Product URL: https://www.faccts.de/orca
Software source components: N/A

Additional Notes:

KIN66 single-point calculations were performed for the complex, isolated
protein fragment, and isolated ligand at identical coordinates. PLA15
full-system DFT calculations used the same three levels of theory at the
original PLA15 geometries. The ORCA 6 method lines were:

B97-3c

    ! b97-3c tightscf slowconv notrah

wB97M-V/def2-TZVPD

    ! wb97m-V def2-TZVPD SCNL tightscf slowconv notrah

wB97M-D3(BJ)/def2-TZVPP

    ! wB97M-D3BJ def2-TZVPP tightscf slowconv notrah


2. Software-specific information:

Name: OpenMM
Version: 8.1.1
System Requirements: CUDA-capable GPU was used in this study
Open Source? (Y/N): Y
Executable URL: N/A
Source Repository URL: https://github.com/openmm/openmm
Developer: OpenMM contributors
Product URL: https://openmm.org
Software source components: N/A

Additional Notes:

Except for datasets based directly on experimentally resolved co-crystal
poses, protein-ligand complexes were subjected to a short restrained
molecular-mechanics minimization of up to 1,000 iterations. Protein and ligand
parameters were ff14SB and GAFF2, respectively. OpenMM 8.1.1 was used with the
OBC2 implicit-solvent model (igb=5 in Amber). Protein-backbone atoms were
positionally restrained; protein side chains and the ligand were allowed to
relax.

After minimization, ligand-centered protein substructures were constructed
using

    R_sphere = R_cut + D_mid

where R_cut = 10 Angstrom and D_mid is the half-diagonal of the ligand bounding
box. Residues were selected as complete units. Neighboring residues and short
sequence gaps were retained to preserve local chain continuity, and open
backbone valences at truncation boundaries were capped with hydrogen atoms.

HiQBind is the co-crystal exception to the restrained-minimization protocol.
Its experimentally resolved ligand poses were retained without the short
OpenMM minimization described above. The same ligand-centered spherical
selection and hydrogen-capping procedure was used to construct the released
protein fragments.


# Additional Notes (such as, will this software not run on certain operating
# systems?):

Additional Notes:

The CSV files are plain UTF-8 text and can be read without ORCA, OpenMM, or the
model software. PDB, SDF, and XYZ files are standard molecular-coordinate
formats. Specialized software is needed only to reproduce the calculations.


#
# Equipment: If specialized equipment generated your data, provide equipment
# name, manufacturer, model, and calibration information. Be sure to include
# specialized file-format information in the data dictionary.
#

3. Equipment-specific information:

Manufacturer: NVIDIA
Model: A100 GPU
Embedded Software / Firmware Name: CUDA
Embedded Software / Firmware Version: N/A
Additional Notes: Used for the short OpenMM minimizations. N/A for calibration.


#
# Dates of Data Collection: List the dates and/or times of data collection.
#

4. Date of data collection (single date, range, approximate date; suggested
format YYYYMMDD): N/A
