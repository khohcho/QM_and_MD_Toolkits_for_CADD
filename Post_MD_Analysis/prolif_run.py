from pathlib import Path
from datetime import datetime
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # NEW: for 2D heatmap
import numpy as np               # NEW: for minor calculations

import MDAnalysis as mda
from MDAnalysis.topology import guessers
import pandas as pd
import prolif as plf
from prolif.plotting import barcode

# LigNetwork is optional
try:
    from prolif.plotting.network import LigNetwork
except Exception:
    LigNetwork = None

# ========= Settings =========
TOP = "md.tpr"
TRJ = "md_no_pbc.xtc"
LIG_RESNAME = "LIG"

AROUND_DYNAMIC = 8.0
VICINITY      = 7.0

RUN_STRIDE = 1
MAX_FRAMES = 10001
N_JOBS = 6

USE_CHUNK_FALLBACK = True
CHUNK = 2450

THRESHOLD = 0.15  
# A pragmatic base threshold for 2D network drawing (reduces clutter)
NETWORK_MIN_THR = 0.15 

OUTDIR = Path("outputs"); OUTDIR.mkdir(exist_ok=True)

# ==== vdW radii (including Cl) ====
VDW_RADII = {
    "H": 1.10, "C": 1.70, "N": 1.55, "O": 1.52, "P": 1.80, "S": 1.80,
    "F": 1.47, "Cl": 1.75, "CL": 1.75, "Br": 1.85, "BR": 1.85, "I": 1.98,
    "Na": 2.27, "NA": 2.27, "K": 2.75, "K+": 2.75, "Ca": 2.31, "CA": 2.31,
    "Mg": 1.73, "MG": 1.73, "Zn": 1.39, "ZN": 1.39, "Fe": 1.32, "FE": 1.32,
}

# ========= Interaction set WITHOUT VdW =========
try:
    import prolif.interactions as plf_inter
    CANDIDATE_NAMES = [
        "HBDonor", "HBAcceptor", "MetalAcceptor", "MetalDonor",
        "Hydrophobic", "PiStacking", "Cationic", "Anionic", "CationPi", "PiCation",
    ]
    INTERACTIONS_NO_VDW = [name for name in CANDIDATE_NAMES if hasattr(plf_inter, name)]
    if not INTERACTIONS_NO_VDW:
        INTERACTIONS_NO_VDW = ["HydrogenBond", "Hydrophobic", "PiStacking"]
except Exception:
    INTERACTIONS_NO_VDW = ["HydrogenBond", "Hydrophobic", "PiStacking"]

print("[i] Utilized interactions (excluding VdW):", INTERACTIONS_NO_VDW)

# ========= Helpers =========
def thr_tag(thr: float) -> str:
    return f"{thr*100:.2f}".replace(".", "_") + "pct"

def safe_save_csv(df, outpath: Path):
    try:
        df.to_csv(outpath); print(f"[✓] Written: {outpath}")
    except PermissionError:
        ts = datetime.now().strftime("%Y%m%d_%H%M%S")
        alt = outpath.with_name(outpath.stem + f"_{ts}.csv")
        df.to_csv(alt); print(f"[!] {outpath.name} is locked. {alt.name} written.")

def safe_save_png(fig, outpath: Path):
    try:
        fig.savefig(outpath, dpi=300, bbox_inches="tight"); print(f"[✓] Written: {outpath}")
    except PermissionError:
        ts = datetime.now().strftime("%Y%m%d_%H%M%S")
        alt = outpath.with_name(outpath.stem + f"_{ts}.png")
        fig.savefig(alt, dpi=300, bbox_inches="tight"); print(f"[!] {outpath.name} is locked. {alt.name} written.")

def write_fallback_edges(fp, threshold: float):
    df = fp.to_dataframe()
    occ = (df.mean() * 100.0).rename("percent").reset_index()
    kept = occ[occ["percent"] >= threshold*100.0].sort_values("percent", ascending=False)
    outp = OUTDIR / f"ligand_network_edges_fallback_thr{thr_tag(threshold)}.csv"
    safe_save_csv(kept, outp)

# ---- Topology repair ----
def ensure_topology(u: mda.Universe):
    try:
        u.trajectory[0]
    except Exception:
        pass
    # elements
    try:
        _ = u.atoms.elements
    except Exception:
        elems = guessers.guess_types(u.atoms.names)
        u.add_TopologyAttr("elements", elems)
    # types
    try:
        _ = u.atoms.types
    except Exception:
        u.add_TopologyAttr("types", [str(el) for el in u.atoms.elements])
        
    # --- SAFETY BELT: If TOP extension is .tpr, SKIP bond guessing! ---
    is_tpr_used = TOP.endswith('.tpr') if 'TOP' in globals() else False

    # bonds
    need_bonds = False
    try:
        nb = u.atoms.bonds
        if nb == 0:
            need_bonds = True
    except Exception:
        need_bonds = True
        
    if need_bonds and not is_tpr_used:
        bonds = guessers.guess_bonds(u.atoms, u.atoms.positions, vdwradii=VDW_RADII)
        u.add_TopologyAttr("bonds", bonds)
        print(f"[i] Bonds reguessed for topology.")
    elif is_tpr_used:
        print(f"[i] 'guess_bonds' skipped because topology source is .tpr (Already embedded).")

    print(f"[i] Topology repair -> elements=OK, types=OK, bonds={u.atoms.bonds}")
    
# --- PDB CLEANER ---
def create_clean_pdb(orijinal_pdb: str, temiz_pdb: str):
    if not os.path.exists(orijinal_pdb):
        print(f"[X] ERROR: Original PDB ({orijinal_pdb}) not found! Please check the filename.")
        return False
        
    print(f"[i] Filtering ATOM lines in {orijinal_pdb} to create {temiz_pdb}...")
    
    with open(orijinal_pdb, "r") as f_in, open(temiz_pdb, "w") as f_out:
        for line in f_in:
            if line.startswith("ATOM"):
                f_out.write(line)
                
    print("[✓] PDB cleaning successful!")
    return True
   
# --- FLAWLESS MAPPER (Works based on Interactive Choice) ---
def map_pdb_to_tpr(pdb_path: str, u_tpr, has_nucleic: bool):
    print(f"[i] Reading reference PDB: {pdb_path}")
    u_pdb = mda.Universe(pdb_path)
    
    # FIRST CRITICAL POINT: The pool is shaped based on the user's 'y' or 'n' answer
    sel_query = "protein or nucleic" if has_nucleic else "protein"
    
    pdb_res = u_pdb.select_atoms(sel_query).residues
    tpr_res = u_tpr.select_atoms(sel_query).residues
    
    pdb_idx = 0
    mapped_count = 0
    
    for t_res in tpr_res:
        if t_res.resname in ["ACE", "NME"]:
            continue
            
        matched = False
        while pdb_idx < len(pdb_res):
            p_res = pdb_res[pdb_idx]
            if t_res.resname == p_res.resname:
                t_res.resid = p_res.resid
                try:
                    t_res.chainID = p_res.chainID
                except AttributeError:
                    pass
                
                pdb_idx += 1
                mapped_count += 1
                matched = True
                break
            else:
                pdb_idx += 1
                
        if not matched:
            print(f"[!] WARNING: {t_res.resname}{t_res.resid} in TPR not found in PDB!")
            
    print(f"[✓] Mapping Complete: {mapped_count}/{len(tpr_res)} residues matched!")   

def run_in_chunks(u, ligand, prot_dyn, stop, vicinity):
    print("[i] Chunk fallback active.")
    all_ifp = {}
    fp_main = plf.Fingerprint(interactions=INTERACTIONS_NO_VDW, vicinity_cutoff=vicinity)
    first = True
    for start in range(0, stop, CHUNK):
        end = min(start + CHUNK, stop)
        print(f"[i] Chunk {start}:{end} (size={end-start})  vicinity={vicinity} Å")
        fp_chunk = plf.Fingerprint(interactions=INTERACTIONS_NO_VDW, vicinity_cutoff=vicinity)
        fp_chunk.run(u.trajectory[start:end:RUN_STRIDE], ligand, prot_dyn, n_jobs=1)
        if first:
            all_ifp.update(fp_chunk.ifp); fp_main = fp_chunk; first = False
        else:
            all_ifp.update(fp_chunk.ifp)
    fp_main.ifp = all_ifp
    return fp_main

def run_once(u, ligand, prot_dyn, stop, vicinity, n_jobs):
    fp = plf.Fingerprint(interactions=INTERACTIONS_NO_VDW, vicinity_cutoff=vicinity)
    frames_slice = u.trajectory[0:stop:RUN_STRIDE]
    print(f"[i] ProLIF run (no VdW): vicinity={vicinity} Å, frames 0:{stop}:{RUN_STRIDE}, n_jobs={n_jobs}")
    try:
        fp.run(frames_slice, ligand, prot_dyn, n_jobs=n_jobs)
    except Exception as e:
        print(f"[!] Error in parallel/serial attempt: {e}\n[i] Trying Serial n_jobs=1...")
        try:
            fp.run(frames_slice, ligand, prot_dyn, n_jobs=1)
        except Exception as e2:
            if USE_CHUNK_FALLBACK:
                print(f"[!] Error in Serial as well ({e2}). Chunk fallback...")
                fp = run_in_chunks(u, ligand, prot_dyn, stop, vicinity)
            else:
                raise
    return fp

# ========= Main Flow =========
def main():
    # --- INTERACTIVE QUERY RUNNING AT THE VERY BEGINNING ---
    print(f"[i] Files: {TOP}, {TRJ}")
    user_input = input("Is there nucleic acid (DNA/RNA) in the system? [y/n]: ").strip().lower()
    has_nucleic = user_input in ['y', 'yes', 'evet']
    
    if has_nucleic:
        print("[i] Nucleic acid mode ACTIVE (Protein + Nucleic will be scanned).")
    else:
        print("[i] Nucleic acid mode CLOSED (Only Protein will be scanned).")

    u = mda.Universe(TOP, TRJ)
    
    # --- Guess missing atom elements for RDKit ---
    for atom in u.atoms:
        if not atom.element:
            atom.element = guessers.guess_atom_element(atom.name)

    # --- FIX: chainID / segid normalization ---
    try:
        ch = list(u.atoms.chainIDs)
    except Exception:
        ch = None
    if ch is None:
        try:
            u.add_TopologyAttr("chainIDs", ["A"] * u.atoms.n_atoms)
        except Exception:
            pass
    else:
        new_ch = [("A" if (c is None or (isinstance(c, str) and c.strip() == "")) else str(c)) for c in ch]
        try:
            u.atoms.chainIDs = new_ch
        except Exception:
            pass
    try:
        sg = list(u.atoms.segids)
    except Exception:
        sg = None
    if sg is None:
        try:
            u.add_TopologyAttr("segids", ["SEG"] * u.atoms.n_atoms)
        except Exception:
            pass
    else:
        new_sg = [("SEG" if (s is None or (isinstance(s, str) and s.strip() == "")) else str(s)) for s in sg]
        try:
            u.atoms.segids = new_sg
        except Exception:
            pass

    # Topology repair
    ensure_topology(u)
    
    # --- PDB CLEANING AND MAPPING PATCH ---
    ORIJINAL_PDB = "PRO.pdb"
    TEMIZ_PDB = "PRO_prolif.pdb"

    if create_clean_pdb(ORIJINAL_PDB, TEMIZ_PDB):
        # SECOND CRITICAL POINT: Preference is passed here and the mapper is triggered
        map_pdb_to_tpr(TEMIZ_PDB, u, has_nucleic)
    else:
        print("[!] WARNING: Mapping skipped because PDB could not be cleaned! Indices will remain in .tpr format.")

    ligand = u.select_atoms(f"resname {LIG_RESNAME}")
    if ligand.n_atoms == 0:
        raise SystemExit(f"[x] Ligand not found: resname {LIG_RESNAME}. "
                         f"Possible names: {sorted(set(r.resname for r in u.residues))[:12]}")

    total = len(u.trajectory)
    stop = min(total, MAX_FRAMES)
    print(f"[i] Frames total={total} | using first {stop} frames")

    # =========================================================================
    # --- DYNAMIC POCKET (PROTEIN + HETATM / ION / WATER) SELECTION ---
    # =========================================================================
    # INFORMATION AND GUIDE (In case you want to add things later):
    # - To search for Zinc: "or resname ZN" or "or resname ZINC"
    # - To search for Iron / HEME: "or resname HEM" or "or resname FE"
    # - To search for Water: "or resname HOH SOL WAT"
    # Example usage: prot_dyn = u.select_atoms(f"(protein or nucleic or resname ZN or resname HOH) and byres around {AROUND_DYNAMIC:.1f} group ligand", ligand=ligand, updating=True)
    # =========================================================================
    
    # THIRD CRITICAL POINT: if 'y', base selection is "protein or nucleic", if 'n' it is only "protein"
    base_prot_sel = "protein or nucleic" if has_nucleic else "protein"
    
    prot_dyn = u.select_atoms(
        f"({base_prot_sel}) and not resname {LIG_RESNAME} and byres around {AROUND_DYNAMIC:.1f} group ligand",
        ligand=ligand, updating=True
    )

    # H count (diagnostic)
    u.trajectory[0]
    nH_lig = ligand.select_atoms("name H*").n_atoms
    nH_prot = u.select_atoms(base_prot_sel + " and name H*").n_atoms
    print(f"[i] H atoms -> ligand: {nH_lig}, protein/nucleic: {nH_prot}")

    # Single run (no VdW)
    fp = run_once(u, ligand, prot_dyn, stop, VICINITY, N_JOBS)

    # ---- Outputs ----
    df = fp.to_dataframe()
    safe_save_csv(df, OUTDIR / "ifp_by_frame.csv")

    keep_cols = df.columns[(df.mean() >= THRESHOLD)]
    df_thr = df.loc[:, keep_cols]
    safe_save_csv(df_thr, OUTDIR / f"ifp_by_frame_thr{thr_tag(THRESHOLD)}.csv")
    print(f"[i] ifp_by_frame: full={df.shape}, thr{thr_tag(THRESHOLD)}={df_thr.shape}")

    occ_all = (df.mean().sort_values(ascending=False).to_frame(name="%") * 100)
    safe_save_csv(occ_all, OUTDIR / "ifp_occurrence_percent_all.csv")
    occ_all_thr = occ_all[occ_all["%"] >= THRESHOLD*100]
    safe_save_csv(occ_all_thr, OUTDIR / f"ifp_occurrence_percent_all_thr{thr_tag(THRESHOLD)}.csv")

    occ_res = (
        df.T.groupby(level=["ligand","protein"]).sum()
          .T.astype(bool).mean().sort_values(ascending=False).to_frame(name="%") * 100
    )
    safe_save_csv(occ_res, OUTDIR / "ifp_occurrence_percent_per_residue.csv")
    occ_res_thr = occ_res[occ_res["%"] >= THRESHOLD*100]
    safe_save_csv(occ_res_thr, OUTDIR / f"ifp_occurrence_percent_per_residue_thr{thr_tag(THRESHOLD)}.csv")

    # Barcode (time axis)
    try:
        ax = barcode.Barcode.from_fingerprint(fp).display(figsize=(12,8))
        times_ns = []
        for ts in u.trajectory[0:stop:RUN_STRIDE]:
            times_ns.append(ts.time * 1e-3)

        n = len(times_ns)
        ax.set_xlabel("Time (ns)")

        num_ticks = 6
        tick_idx = np.linspace(0, n - 1, num_ticks, dtype=int)
        ax.set_xticks(tick_idx)
        ax.set_xticklabels([f"{times_ns[i]:.1f}" for i in tick_idx])
        
        safe_save_png(ax.figure, OUTDIR / "ifp_barcode.png")
    except Exception as e:
        print(f"[!] Barcode could not be drawn: {e}")

    # ==== 2D interaction network (HTML) ====
    try:
        lig_mol = plf.Molecule.from_mda(ligand)
        TH2D = max(THRESHOLD, NETWORK_MIN_THR)  # reduce visual clutter
        try:
            view = fp.plot_lignetwork(lig_mol, threshold=TH2D)
            view.save(str(OUTDIR / f"lignetwork_thr{thr_tag(TH2D)}.html"))
            print(f"[✓] Written: {OUTDIR/f'lignetwork_thr{thr_tag(TH2D)}.html'}")
        except Exception as e:
            print(f"[i] fp.plot_lignetwork could not be used ({e}) → LigNetwork fallback.")
            if LigNetwork is not None:
                net = LigNetwork.from_fingerprint(fp, lig_mol, kind="aggregate", threshold=TH2D)
                net.save(str(OUTDIR / f"lignetwork_thr{thr_tag(TH2D)}.html"))
                print(f"[✓] Written: {OUTDIR/f'lignetwork_thr{thr_tag(TH2D)}.html'}")
            else:
                print("[!] No LigNetwork; CSV fallback written.")
                write_fallback_edges(fp, THRESHOLD)

        try:
            top_frames = df.sum(axis=1).sort_values(ascending=False).head(3).index.tolist()
            if LigNetwork is not None and len(top_frames) > 0:
                for f in top_frames:
                    net_f = LigNetwork.from_fingerprint(fp, lig_mol, kind="frame", frame=f, threshold=0.0)
                    net_f.save(str(OUTDIR / f"lignetwork_frame{f}.html"))
                    print(f"[✓] Written: {OUTDIR/f'lignetwork_frame{f}.html'}")
        except Exception as e:
            print(f"[!] Frame-based 2D network could not be generated: {e}")

    except Exception as e:
        print(f"[!] 2D network error: {e}")

    # ==== 3D VISUALIZATION (py3Dmol / Complex3D) ====
    try:
        try:
            df_ifp = fp.to_dataframe()
            TOP_FRAMES_N = 3
            frames_to_save = df_ifp.sum(axis=1).sort_values(ascending=False).head(TOP_FRAMES_N).index.tolist()
        except Exception:
            frames_to_save = [0]

        use_segid = getattr(fp, "use_segid", False)

        for f in frames_to_save:
            u.trajectory[f]
            lig_mol  = plf.Molecule.from_mda(ligand,   use_segid=use_segid)
            prot_mol = plf.Molecule.from_mda(prot_dyn, use_segid=use_segid)

            try:
                view3d = fp.plot_3d(lig_mol, prot_mol, frame=f, display_all=False)
            except TypeError:
                view3d = fp.plot_3d(lig_mol, prot_mol, None, frame=f, display_all=False)

            out3d = OUTDIR / f"complex3d_frame{f}.html"
            with open(out3d, "w", encoding="utf-8") as html_dosyasi:
                html_dosyasi.write(view3d._make_html())
            
            print(f"[✓] 3D saved: {out3d.absolute()}")
    except Exception as e:
        print(f"[!] 3D visualization failed: {e}")


if __name__ == "__main__":
    main()
