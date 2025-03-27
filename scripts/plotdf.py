#!/usr/bin/env python3

import os, re, glob, argparse, warnings
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import mplcursors

warnings.simplefilter('ignore')

# ---- Argument Parser ----
parser = argparse.ArgumentParser(
    prog='plotdf.py',
    usage='%(prog)s <lc_file> [--derivatives dF_type1 dF_type2 ...]',
    description='Interactive visualization of microlensing light curves and Fisher matrix derivatives.'
)
parser.add_argument('lc_file', type=str, help="Directory containing light curve files")
parser.add_argument('--derivatives', nargs='+', default=["parameter_derivatives"],
                    help="List of derivative types to include")
args = parser.parse_args()

# ---- Config ----
header = [
    "Simulation_time", "measured_relative_flux", "measured_relative_flux_error", "true_relative_flux",
    "true_relative_flux_error", "observatory_code", "saturation_flag", "best_single_lens_fit", 
    "parallax_shift_t", "parallax_shift_u", "BJD", "source_x", "source_y", "lens1_x", "lens1_y", 
    "lens2_x", "lens2_y", "x", "y", "z"
]
parameters = ["t0", "tE", "u0", "alpha", "s", "q", "rs", "piEN", "piEE", "Fbase0", "fs0", "Fbase1", "fs1", "Fbase2", "fs2"]
outparameters = ["t0lens1", "tE_ref", "u0lens1", "alpha", "Planet_s", "Planet_q", "rho", "piEN", "piEE", "Obs_0_fs", "Obs_1_fs", "Obs_2_fs"]
log_parameters = ["tE", "s", "q", "rs"]
for dtype in args.derivatives:
    header.extend([f"{dtype}_{p}" for p in parameters])

lc_files = glob.glob(os.path.join(args.lc_file))
if not lc_files:
    raise SystemExit("No .lc files found. Exiting.")

# Derive .out file from the first .lc file
match = re.match(r'^(.*)_(\d+)\.det\.lc$', os.path.basename(lc_files[0]))
if not match:
    raise ValueError(f"Invalid lc file format: {os.path.basename(lc_files[0])}")

out_file = os.path.join(os.path.dirname(os.path.basename(lc_files[0])), match.group(1) + ".out")
fisher_data = pd.read_csv(out_file, delim_whitespace=True)
if not out_file:
	raise SystemExit("No .out file found. Exiting.")

event_id = int(match.group(2))
fisher_error_dict = {}

# ---- Fisher Matrix Construction ----
def compute_fisher(df, dtype):
    cols = [c for c in df.columns if c.startswith(dtype)]
    if not cols: return None
    err2 = df['measured_relative_flux_error'] ** 2
    fm = np.array([[np.sum(df[c1]*df[c2]/err2) for c2 in cols] for c1 in cols])
    try:
        cov = np.linalg.pinv(fm)
        return dict(zip(cols, np.sqrt(np.diag(cov))))
    except np.linalg.LinAlgError:
        return dict(zip(cols, ["N/A"] * len(cols)))

for lc in lc_files:
    try:
        df = pd.read_csv(lc, delim_whitespace=True, comment='#', header=None, names=header)
        for dtype in args.derivatives:
            res = compute_fisher(df, dtype)
            if res: fisher_error_dict.update(res)
    except FileNotFoundError:
        print(f"File not found: {lc}, skipping.")

# ---- Print Table ----
if event_id and 'EventID' in fisher_data.columns:
    row = fisher_data[fisher_data['EventID'] == event_id].squeeze()
    if not row.empty:
        print(f"\nEventID = {event_id}")
        header_row = "Parameter         | Value         | OutFile_Unc         | " + " | ".join(f"FM_{d}" for d in args.derivatives)
        print(header_row + "\n" + "-" * len(header_row))
        for out_p, p in zip(outparameters, parameters):
            val = row.get(out_p, "N/A")
            unc = row.get(f"ObsGroup_0_{p}_err", "N/A")
            if p in log_parameters and isinstance(val, (int, float)) and isinstance(unc, (int, float)):
                unc *= val * np.log(10)
            ests = []
            for dtype in args.derivatives:
                k = f"{dtype}_{p}"
                fm = fisher_error_dict.get(k, "N/A")
                if p in log_parameters and isinstance(fm, (int, float)) and isinstance(val, (int, float)):
                    fm *= val * np.log(10)
                ests.append(f"{fm:<10.6f}" if isinstance(fm, (int, float)) else f"{fm:<10}")
            val_str = f"{val:<13.6f}" if isinstance(val, (int, float)) else f"{val:<13}"
            unc_str = f"{unc:<18.6f}" if isinstance(unc, (int, float)) else f"{unc:<18}"
            print(f"{p:<17} | {val_str} | {unc_str} | " + " | ".join(ests))
    else:
        print(f"EventID {event_id} not found.")

# ---- Plot Lightcurve and Derivatives ----
def plot_lc_derivatives(lc):
    try:
        df = pd.read_csv(lc, delim_whitespace=True, comment='#', header=None, names=header)
        df = df[df['observatory_code'] == 0]
        fsm = fsm = pd.read_csv(lc, delim_whitespace=True, header=None, comment=None, engine='python', nrows=4, index_col=False)
        fs, m0 = fsm.iloc[0, 1], fsm.iloc[-1, 1] + 2.5 * np.log10(fsm.iloc[0, 1])
        mi = m0 - 2.5 * np.log10(fs * df["measured_relative_flux"] + 1 - fs)

        palette = sns.color_palette("tab10", n_colors=len(args.derivatives))
        color_map = {d: palette[i] for i, d in enumerate(args.derivatives)}

        nrows = (len(parameters) + 1) // 2
        fig, axes = plt.subplots(nrows + 1, 2, figsize=(12, (nrows + 1) * 3), sharex=True)
        axes = axes.flatten()
        fig.suptitle(f"Fisher derivatives & light curve - {os.path.basename(lc)}")

        # Light curve
        axes[0].plot(df["Simulation_time"], mi, 'o', markersize=2, color='black', alpha=0.7)
        axes[0].set_ylabel("Magnitude"); axes[0].invert_yaxis(); axes[0].grid()

        # Derivatives
        for i, p in enumerate(parameters, 1):
            for d in args.derivatives:
                col = f"{d}_{p}"
                if col in df.columns:
                    axes[i].plot(df["Simulation_time"], df[col], 'o--', markersize=4,
                                 color=color_map[d], alpha=0.2, label=d)
            axes[i].set_ylabel(p); axes[i].grid()

        # Legend and metadata
        legend_handles = [
    plt.Line2D([0], [0], marker='o', color='black', label="Light Curve", markersize=6)
] + [plt.Line2D([0], [0], marker='o', color='w', markerfacecolor=color, markersize=8, label=d)
     for d, color in color_map.items()]
        
        fig.legend(handles=legend_handles, loc='lower left', ncol=len(color_map), fontsize='small')

        if event_id and 'EventID' in fisher_data.columns:
            row = fisher_data[fisher_data['EventID'] == event_id].squeeze()
            if not row.empty:
                param_text = " | ".join(
                    f"{p:<6}: {row.get(p, 'N/A'):.4f}" if isinstance(row.get(p), (float, int)) else f"{p:<6}: N/A"
                    for p in outparameters)
                fig.text(0.5, 0.01, f"Fisher Parameters:\n{param_text}", ha='center', fontsize='large')

        
        # Hide unused axes
        for j in range(len(parameters) + 1, len(axes)):
            fig.delaxes(axes[j])
            
        
        valid_axes = [ax for ax in axes if ax in fig.axes]# second to last valid subplot
        bottom_ax = valid_axes[-1]   # Last valid subplot
        second_bottom_ax = valid_axes[-2]  # Second to last valid subplot
        bottom_ax.tick_params(labelbottom=True)
        second_bottom_ax.tick_params(labelbottom=True)
        second_bottom_ax.set_xlabel("Simulation Time")
        bottom_ax.set_xlabel("Simulation Time")

        mplcursors.cursor(axes, hover=True).connect("add", lambda sel:
            sel.annotation.set_text(f"Time: {sel.target[0]:.8f}\nValue: {sel.target[1]:.8f}"))

        plt.tight_layout(rect=[0, 0.02, 1, 0.98])
        plt.show()

    except FileNotFoundError:
        print(f"File not found: {lc}, skipping.")

for lc in lc_files:
    plot_lc_derivatives(lc)
