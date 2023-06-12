# DockingProtocol

Repository to create a concrete workflow to work with in drug discovery projects.


## Workflow

Target
            --------->             --------->  Docking  --------->  tPELE  ---------> Analysis
Inhibitors               LigPrep
    
### Input

1. Target to inhibit prepared and in format pdb, mol2 or sdf.
2. Csv file with only SMILES and id of all the inhibitor ligands.

### Ligprep

Generation of tautomers and isomers as well as protonation of the inhibitor ligands in selected pH +- pH_tolerance.

### Docking

Docking of the inhibitors to the target.

### tPELE

Refinement of the docking pose obtained in the previous step in the pipeline.

### Analysis

Part to analyze the results obteined in the tPELE simulations.
