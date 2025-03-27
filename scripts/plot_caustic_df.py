import os, re, glob, argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import mplcursors
import VBMicrolensing


VBM = VBMicrolensing.VBMicrolensing()
parser = argparse.ArgumentParser(
    prog='plot_caustic_df.py',
    usage='%(prog)s <lc_file>',
    description="Visualize caustics"
)
parser.add_argument("lcfile", type=str, help="lightcurve files (*.lc)")
args = parser.parse_args()

lc_files = glob.glob(args.lcfile)
if not lc_files:
    raise SystemExit("No .lc files found. Exiting.")

header = [
    "Simulation_time", "measured_relative_flux", "measured_relative_flux_error", "true_relative_flux",
    "true_relative_flux_error", "observatory_code", "saturation_flag", "best_single_lens_fit",
    "parallax_shift_t", "parallax_shift_u", "BJD", "source_x", "source_y",
    "lens1_x", "lens1_y", "lens2_x", "lens2_y", "x", "y", "z"
]

parameters = ["t0", "tE", "u0", "alpha", "s", "q", "rs", "piEN", "piEE", "Fbase0", "fs0", "Fbase1", "fs1", "Fbase2", "fs2"]
outparams = ["t0lens1", "tE_ref", "u0lens1", "alpha", "Planet_s", "Planet_q", "rho", "piEN", "piEE", "Obs_0_fs", "Obs_1_fs", "Obs_2_fs"]
log_params = {"tE", "s", "q", "rs"}
derivative_types = ["dF_delta4", "dF_delta2", "dF_delta10"]
header += [f"{d}_{p}" for d in derivative_types for p in parameters]

# --- Get EventID and Fisher Row ---
filename = os.path.basename(lc_files[0])
match = re.match(r'^(.*)_(\d+)\.det\.lc$', os.path.basename(lc_files[0]))
if not match:
    raise ValueError(f"Invalid lc file format: {os.path.basename(lc_files[0])}")

out_file = os.path.join(os.path.dirname(os.path.basename(lc_files[0])), match.group(1) + ".out")
fisher_data = pd.read_csv(out_file, delim_whitespace=True)
if not out_file:
	raise SystemExit("No .out file found. Exiting.")

event_id = int(match.group(2))

fisher_row = fisher_data.loc[fisher_data["EventID"] == event_id].squeeze() if event_id else None
if fisher_row is None or fisher_row.empty:
    raise ValueError(f"EventID {event_id} not found in {out_file}")

# --- Print parameter values and uncertainties ---
print(f"EventID = {event_id}\n{'Parameter':10} | {'Value':>10} | {'FM_Uncertainty'}\n" + "-" * 35)
for p_out, p in zip(outparams, parameters):
    val = fisher_row.get(p_out, "N/A")
    err = fisher_row.get(f"ObsGroup_0_{p}_err", "N/A")
    if p in log_params and isinstance(val, (int, float)) and isinstance(err, (int, float)):
        err *= val * np.log(10)
    print(f"{p:10} | {val:10} | {err}")

# --- Loop through LC files ---
for lc in lc_files:
    try:
        df = pd.read_csv(lc, delim_whitespace=True, comment="#", header=None, names=header)
        time = df["Simulation_time"].values
        y1, y2 = df["source_x"].values, df["source_y"].values

        # Prepare binary microlensing parameters
        t0, tE = fisher_row["t0lens1"], fisher_row["tE_ref"]
        s, q = fisher_row["Planet_s"], fisher_row["Planet_q"]
        params = np.log([s, q, fisher_row["rho"], tE]).tolist()
        params.insert(2, fisher_row["u0lens1"])
        params.insert(3, fisher_row["alpha"])
        params.append(t0)

        # Compute light curve + caustics
        mag = VBM.BinaryLightCurve(params, time)
        caustics = VBM.Caustics(s, q)

        # Plot
        fig, ax = plt.subplots(figsize=(6, 6))
        all_x, all_y = np.concatenate([c[0] for c in caustics]), np.concatenate([c[1] for c in caustics])
        for cx, cy in caustics:
            ax.plot(cx, cy, 'k', lw=3)

        # Source trajectory
        sc = ax.scatter(mag[1], mag[2], c=time, cmap='coolwarm', s=1)
        plt.colorbar(sc, ax=ax, label="Time (days)")
        ax.set(xlabel="x", ylabel="y", title="Caustic and Source Trajectory", aspect='equal')
        cursor = mplcursors.cursor(sc, hover=True)
        cursor.connect("add", lambda sel: sel.annotation.set_text(f"Time: {time[sel.index]:.4f}"))


        mx, my = np.mean(mag[1]), np.mean(mag[2])
        dx, dy = mag[1][-1] - mag[1][0], mag[2][-1] - mag[2][0]
        ax.annotate('', xy=(mx + dx * 0.002, my + dy * 0.002), xytext=(mx, my),
                    arrowprops=dict(arrowstyle='->', color='red', linewidth=2.5))
        
        #plt.ylim(-1.0, 1.0)
        ax.legend()
        plt.show()

    except FileNotFoundError:
        print(f"File not found: {lc}, skipping.")
