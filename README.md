# DockingProtocol 🧬

A Python-based workflow to streamline **ligand preparation, docking, refinement, and analysis** in **drug discovery projects**.  
This repository integrates ligand preparation, docking (Glide, rDock, EquiBind), and PELE refinement to create a consistent and automated virtual screening pipeline.

---

## 🚀 Overview

`DockingProtocol` helps computational chemists and structural biologists build reproducible docking workflows.  
It covers all stages — from ligand preparation to post-simulation analysis — using industry-standard and research-grade tools.

### Workflow Diagram

```mermaid
graph LR
A[Inhibitors] --> B[LigPrep]
B --> C[Docking]
D[Target] --> C[Docking]
C -- Best poses --> E[PELE]
E --> F[Analysis]
