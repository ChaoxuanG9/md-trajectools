# residence_time.py
"""
Author: Chaoxuan Gu (chaoxuan_gu@brown.com)
GitHub: https://github.com/ChaoxuanG9/md-trajectools
License: MIT

If you use this script in your work, please cite the repo.
"""
# ----------------------------------------------
# General-purpose Residence Time Analysis Script
# Supports custom .xyz format trajectory files
# Works in two region modes: EDL (slab) and Shell (solvation)
# Allows atom selection by charge, label, or name
# Input via CLI or YAML config file
# ----------------------------------------------

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import glob
import os
import argparse
import logging
import yaml
from datetime import datetime

# ----------------------------
# Load YAML configuration
# ----------------------------
def load_config_yaml(path):
    with open(path, 'r') as f:
        return yaml.safe_load(f)

# ----------------------------
# Argument parser with optional YAML support
# CLI > YAML > default
# ----------------------------
def parse_args():
    parser = argparse.ArgumentParser(description="General-purpose Residence Time Analysis")
    parser.add_argument('--config', help="Path to YAML config file with all arguments")
    parser.add_argument('--path', help="Glob path to trajectory .xyz files")
    parser.add_argument('--mode', choices=['edl', 'shell'], default='edl', help="How to define region of interest")
    parser.add_argument('--zmin', type=float, help="Z lower bound (for EDL mode)")
    parser.add_argument('--zmax', type=float, help="Z upper bound (for EDL mode)")
    parser.add_argument('--cutoff', type=float, help="Cutoff distance (for shell mode)")
    parser.add_argument('--exclude', type=int, help="Number of initial frames to skip")
    parser.add_argument('--outdir', help="Output directory for results")
    parser.add_argument('--from0', action='store_true', help="Use R_i(t) · R_i(0) instead of time-averaged autocorrelation")
    parser.add_argument('--select-method', choices=['charge', 'label', 'name'], help="How to identify atoms of interest")
    parser.add_argument('--select-value', help="Value for atom selection (e.g., -0.32 or 'N')")
    parser.add_argument('--box-x', type=float, help="PBC cell length in X")
    parser.add_argument('--box-y', type=float, help="PBC cell length in Y")
    args = parser.parse_args()

    # Load defaults from YAML if provided
    if args.config:
        yaml_config = load_config_yaml(args.config)
        for key, val in yaml_config.items():
            if getattr(args, key, None) is None:
                setattr(args, key, val)
    return args

# ----------------------------
# Logging to file and terminal
# ----------------------------
def setup_logging(outdir):
    os.makedirs(outdir, exist_ok=True)
    log_path = os.path.join(outdir, "residence_time.log")
    logging.basicConfig(
        filename=log_path,
        level=logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s"
    )

# ----------------------------
# Atom selection utility
# ----------------------------
def select_atoms(df, method, value):
    """Select atoms from a frame using label, name, or charge"""
    if method == "charge":
        return df[df['charge'] == float(value)]
    elif method == "label":
        return df[df['label'] == value]
    elif method == "name":
        return df[df['name'] == value]
    else:
        raise ValueError(f"Invalid selection method: {method}")

# ----------------------------
# Trajectory loader from .xyz files
# ----------------------------
def load_trajectories(path_pattern):
    """Load all matching .xyz trajectory files into a dict of DataFrames"""
    trj_df = {}
    for path in glob.glob(path_pattern):
        key = os.path.basename(path).replace(".txt", "")
        logging.info(f"Loading {key} from {path}")
        trj_df[key] = pd.read_csv(path, delim_whitespace=True, comment="#")
    return trj_df

# ----------------------------
# Minimum image distance with PBC in x/y
# ----------------------------
def calc_distance(coord1, coord2, box_x, box_y):
    shifts = [-1, 0, 1]
    d_min = float('inf')
    for dx in shifts:
        for dy in shifts:
            translated = np.array([
                coord2[0] + dx * box_x,
                coord2[1] + dy * box_y,
                coord2[2]
            ])
            dist = np.linalg.norm(np.array(coord1) - translated)
            d_min = min(d_min, dist)
    return d_min

# ----------------------------
# Create a region-checking function
# This returns True if an atom is in the region of interest
# ----------------------------
def make_region_checker(args, trj_df):
    if args.mode == "edl":
        # Simple slab region along z-axis
        def region_func(atom, frame_num=None):
            return args.zmin <= atom['Z'] <= args.zmax
        return region_func

    elif args.mode == "shell":
        # For shell mode, find centers (e.g., N atoms) first
        ion_centers = {}
        for sys_df in trj_df.values():
            for frame in sys_df['frameNum'].unique():
                snap = sys_df[sys_df['frameNum'] == frame]
                ions = snap[snap['label'] == args.select_value]
                ion_centers[frame] = ions[['X', 'Y', 'Z']].values

        # For each frame, check if atom is within cutoff from any ion
        def region_func(atom, frame_num):
            coord = atom[['X', 'Y', 'Z']].values
            for center in ion_centers.get(frame_num, []):
                if calc_distance(coord, center, args.box_x, args.box_y) < args.cutoff:
                    return True
            return False

        return region_func

# ----------------------------
# Find time when correlation decays to 1/e
# ----------------------------
def find_decay_time_to_e_inverse(corr, time):
    target = 1 / np.e
    for i in range(1, len(corr)):
        if corr[i] <= target:
            t0, t1 = time[i - 1], time[i]
            c0, c1 = corr[i - 1], corr[i]
            return t0 + (target - c0) * (t1 - t0) / (c1 - c0)
    return np.nan

# ----------------------------
# Compute autocorrelation of residence binary function
# ----------------------------
def compute_residence_correlation(R_matrix, from0):
    n_frames, n_mol = R_matrix.shape
    corr = []
    if from0:
        ref = R_matrix[0]
        for t in range(n_frames):
            corr.append(np.dot(R_matrix[t], ref) / n_mol)
    else:
        for tau in range(n_frames):
            values = [np.dot(R_matrix[t], R_matrix[t + tau]) / n_mol
                      for t in range(n_frames - tau)]
            corr.append(np.mean(values))
    return np.array(corr)

# ----------------------------
# Save plot of autocorrelation
# ----------------------------
def plot_correlation(time, corr, label, out_path):
    plt.figure(figsize=(5, 4))
    plt.plot(time, corr, label=label)
    plt.xlabel("Time (ps)")
    plt.ylabel("Autocorrelation")
    plt.title(f"Residence Correlation: {label}")
    plt.legend()
    plt.tight_layout()
    plt.savefig(out_path, dpi=300)
    plt.close()

# ----------------------------
# Main analysis pipeline (works for both EDL and Shell)
# ----------------------------
def analyze(trj_df, args):
    os.makedirs(args.outdir, exist_ok=True)
    region_func = make_region_checker(args, trj_df)
    decay_summary = []

    for sysname, df in trj_df.items():
        logging.info(f"Analyzing system: {sysname}")
        frames = sorted(df['frameNum'].unique())
        dt = 5  # timestep in ps
        n_frames = len(frames) - args.exclude
        time = np.arange(0, n_frames * dt, dt)
        R_matrix = []

        for idx in range(args.exclude, len(frames)):
            frame = frames[idx]
            snap = df[df['frameNum'] == frame]
            selected = select_atoms(snap, args.select_method, args.select_value)
            R_vec = [int(region_func(row, frame)) for _, row in selected.iterrows()]
            R_matrix.append(R_vec)

        R_matrix = np.array(R_matrix)
        corr = compute_residence_correlation(R_matrix, args.from0)
        tau = find_decay_time_to_e_inverse(corr, time)

        # Save results
        df_out = pd.DataFrame({'time': time, 'corr': corr})
        df_out.to_csv(os.path.join(args.outdir, f"{sysname}_res_time.csv"), index=False)
        plot_correlation(time, corr, sysname, os.path.join(args.outdir, f"{sysname}_res_time.png"))
        decay_summary.append((sysname, tau))
        logging.info(f"{sysname}: tau(1/e) = {tau:.2f} ps")

    # Save summary of decay times
    pd.DataFrame(decay_summary, columns=["system", "tau_e"]).to_csv(
        os.path.join(args.outdir, "decay_times.csv"), index=False
    )

# ----------------------------
# Entrypoint
# ----------------------------
def main():
    args = parse_args()
    setup_logging(args.outdir)
    logging.info("Starting analysis with args:")
    logging.info(vars(args))
    trj_df = load_trajectories(args.path)
    analyze(trj_df, args)

if __name__ == "__main__":
    main()
