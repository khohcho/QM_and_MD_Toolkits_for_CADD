# Post-MD Analysis Toolkit

## prolif.py

A professional, robust Python toolkit for analyzing Protein-Ligand and Nucleic Acid-Ligand Interaction Fingerprints (IFP) from Molecular Dynamics (MD) trajectories. This script leverages **MDAnalysis** and **ProLIF** to extract, visualize, and quantify key pharmacophoric interactions while eliminating background noise.

## Key Features

*   **Interactive Target Selection:** Dynamically supports both standard Protein-Ligand systems and complex Protein-Nucleic Acid-Ligand environments via a simple prompt.
*   **Targeted Interaction Filtering:** Deliberately excludes ubiquitous Van der Waals (VdW) contacts to maximize the signal-to-noise ratio, focusing strictly on critical lock-and-key interactions (Hydrogen bonds, Pi-Stacking, Hydrophobic enclosures, etc.).
*   **Auto-Topology Repair:** Automatically handles missing elements, chain IDs, and segment IDs to prevent common `.xtc`/`.tpr` parsing errors.
*   **Terminal-Safe 3D Rendering:** Bypasses Jupyter-dependent rendering by "forcibly" extracting the underlying HTML/JS code from `py3Dmol`, safely writing 3D interaction maps to your local disk even on headless Ubuntu terminal sessions.
*   **PDB Mapping Patch:** Automatically cleans and maps your reference `.pdb` files to your `.tpr` index to ensure perfect residue numbering.

---

## Installation & Recommended Environment

This script is highly sensitive to library versions. It has been extensively tested and verified to work flawlessly with the following environment setup.

### 1. Create and Activate the Conda Environment
```bash
conda create -n prolif_env python=3.12 -y
conda activate prolif_env
```

### 2. Install Required Dependencies
```bash
conda install -c conda-forge mdanalysis prolif pandas matplotlib ipython py3dmol -y
```

### 3. Tested & Verified Package Versions
If you encounter any issues in the future, force-installing these specific versions is highly recommended:

*   python: 3.12.13
*   prolif: 2.1.0
*   mdanalysis: 2.10.0
*   pandas: 3.0.2
*   matplotlib: 3.10.9
*   ipython: 9.13.0
*   py3dmol: 2.5.4

## Required Input Files
Before running the script, ensure the following files are present in your working directory:

*   `md.tpr`: Your MD run input file (topology).
*   `md_no_pbc.xtc`: Your MD trajectory file (make sure Periodic Boundary Conditions are removed and the complex is centered).
*   `PRO.pdb`: **[CRITICAL]** The pure protein structure file (usually generated via `pdb2gmx` or manually extracted). The script requires this exact naming (`PRO.pdb`) to perform internal PDB cleaning and residue mapping (`PRO_prolif.pdb`).

## Usage
Simply execute the script in your terminal:

```bash
python prolif.py
```

The script will prompt you:

```
Is there nucleic acid (DNA/RNA) in the system? [y/n]:
```

Answer according to your system, and the automated pipeline will begin processing your trajectory frames. All outputs (CSVs, Barcodes, 2D HTML Networks, and 3D HTML Complex views) will be securely saved in the automatically generated `outputs/` folder.

## Known Issues & Technical Notes
*   **Deprecation Warnings:** You may see `DeprecationWarning` messages in the terminal output regarding MDAnalysis or ProLIF internals. These are expected due to the rapid evolution of these libraries and do not affect the accuracy of your results.
*   **3D HTML Forced Saves:** The terminal output will indicate `[✓] 3D saved: .../complex3d_frameX.html`. This utilizes a custom bypass method to write standard web code directly to the disk, overcoming the native limitation of `py3Dmol` requiring an interactive Jupyter Notebook. You can open these `.html` files in any standard web browser.

## License
This project is licensed under the GPL-3.0 License to maintain compatibility with ProLIF and support the open-source computational chemistry community.
