
# DockingProtocol

Repository to create a concrete workflow to work with in drug discovery projects.


## Workflow

```mermaid
graph LR
A[Inhibitors] --> B[LigPrep]
B --> C[Docking]
D[Target] --> C[Docking]
C -- Best poses --> E[PELE]
E --> F[Analysis]
```
---

### 1. Input

1. Target to inhibit prepared and in format pdb or mol2.
2. Csv file with only SMILES and id of all the inhibitor ligands.

**- Module:** inputPrepare.py

**- Class:** InputPreparation

**- Methods:** setUpLigPrepJob

#### 1.1. Ligprep

Generation of tautomers and isomers as well as protonation of the inhibitor ligands in selected pH +- pH_tolerance.

---

### 2. Docking

Docking of the inhibitors to the target or rescoring docked poses.

**Modules:** 
1. dockingJob.py 

	**- Class:** DockingJob
	
	**- Methods:** setGlideDocking, setRdockDockingset, setEquibindDocking, rDockRescore, and GlideRescore
	
2. dockingAnalysis.py

	**- Class:** DockingAnalyzer

	**- Methods:** glideAnalysis, rdockAnalysis

---

### 3. PELE

Refinement of the docking pose obtained in the previous step in the pipeline.

**- Module:** peleJob.py

**- Class:** PELE

**- Methods:** setGlideToPELESimulation, setRdockToPELESimulation, setEquibindToPELESimulation, and PELEDownloader.

---

### 4. Analysis

Part to analyze the results obtained in the PELE simulations.

**- Module:** peleAnalysis.py

**- Class:** PELEAnalyzer

**- Methods:** experimentalDataCollector, equibindDataTrimming, PELEDataCollector, and correlationPlotter.

---
