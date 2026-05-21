import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import glob
import re
from pathlib import Path
from scipy.optimize import curve_fit

# ============================================================
# SETTINGS
# ============================================================
file_pattern = "C:/Users/firew/Documents/GitHub/imphys/Analogue/Data/*.csv"   # folder containing CSV files

# Set ONE of these approaches:
use_named_columns = False     # True if your CSV has named columns
time_col = 0
input_col = 1
output_col = 2

# If use_named_columns = True, set these:
time_name = "Time"
input_name = "Input"
output_name = "Output"

# ============================================================
# HELPERS
# ============================================================
def extract_frequency(filename):
    """
    Extract frequency from filename, e.g.:
    '0.5Hz.csv' -> 0.5
    'trace_2Hz.csv' -> 2.0
    """
    match = re.search(r'(\d+\.?\d*)\s*Hz', Path(filename).stem, re.IGNORECASE)
    if match:
        return float(match.group(1))
    raise ValueError(f"Could not extract frequency from filename: {filename}")

def sine_model(t, A, phi, C, f):
    """
    Sine model with fixed frequency f:
    y = A sin(2π f t + φ) + C
    """
    return A * np.sin(2 * np.pi * f * t + phi) + C

def fit_trace(time, signal, freq):
    """
    Fit a sine wave of fixed frequency to noisy data.
    Returns:
        amplitude, amplitude_uncertainty, phase, phase_uncertainty, offset
    """
    time = np.asarray(time, dtype=float)
    signal = np.asarray(signal, dtype=float)

    mask = np.isfinite(time) & np.isfinite(signal)
    time = time[mask]
    signal = signal[mask]

    if len(time) < 5:
        raise ValueError("Not enough valid data points for fitting")

    A_guess = 0.5 * (np.max(signal) - np.min(signal))
    C_guess = np.mean(signal)

    if A_guess == 0:
        phi_guess = 0.0
    else:
        y0 = signal[0] - C_guess
        phi_guess = np.arcsin(np.clip(y0 / A_guess, -1, 1))

    popt, pcov = curve_fit(
        lambda t, A, phi, C: sine_model(t, A, phi, C, freq),
        time,
        signal,
        p0=[A_guess, phi_guess, C_guess],
        maxfev=20000
    )

    A_fit, phi_fit, C_fit = popt
    perr = np.sqrt(np.diag(pcov))

    A_err = perr[0]
    phi_err = perr[1]

    # Force amplitude positive, absorb sign into phase if needed
    if A_fit < 0:
        A_fit = -A_fit
        phi_fit += np.pi

    # Wrap phase into [-pi, pi]
    phi_fit = np.arctan2(np.sin(phi_fit), np.cos(phi_fit))

    return A_fit, A_err, phi_fit, phi_err, C_fit

def bode_amplitude(f, omega0, gamma, K):
    """
    Magnitude response of driven damped oscillator:
    G(f) = K / sqrt((omega0^2 - omega^2)^2 + (gamma*omega)^2)
    where omega = 2πf
    """
    omega = 2 * np.pi * f
    return K / np.sqrt((omega0**2 - omega**2)**2 + (gamma * omega)**2)

def bode_phase(f, omega0, gamma):
    """
    Phase lag of driven damped oscillator in radians.
    """
    omega = 2 * np.pi * f
    return -np.arctan2(gamma * omega, (omega0**2 - omega**2))

def choose_phase_branch(phi_out, phi_in, freq, omega0_ref, gamma_ref):
    """
    Compute output phase relative to input, then choose the equivalent
    phase branch that is closest to the expected damped-oscillator phase.

    Returns:
        dphi_chosen : wrapped phase difference in radians
    """
    dphi_raw = phi_out - phi_in

    candidates = np.array([
        dphi_raw - 2*np.pi,
        dphi_raw - np.pi,
        dphi_raw,
        dphi_raw + np.pi,
        dphi_raw + 2*np.pi
    ])

    phi_expected = bode_phase(freq, omega0_ref, gamma_ref)
    dphi_chosen = candidates[np.argmin(np.abs(candidates - phi_expected))]
    dphi_chosen = np.arctan2(np.sin(dphi_chosen), np.cos(dphi_chosen))

    return dphi_chosen

# ============================================================
# LOAD AND PROCESS FILES
# ============================================================
files = sorted(glob.glob(file_pattern))
if not files:
    raise FileNotFoundError(f"No CSV files found matching pattern: {file_pattern}")

# Rough reference phase-model parameters for branch selection
freq_list_for_guess = [extract_frequency(f) for f in files]
freq_array_guess = np.array(sorted(freq_list_for_guess))
f_ref_guess = np.median(freq_array_guess)
omega0_ref_guess = 2 * np.pi * f_ref_guess
gamma_ref_guess = omega0_ref_guess / 2

frequencies = []
gains = []
gain_unc = []
phases = []
phase_unc = []

for file in files:
    try:
        freq = extract_frequency(file)
        df = pd.read_csv(file)

        if use_named_columns:
            t = df[time_name].to_numpy()
            vin = df[input_name].to_numpy()
            vout = df[output_name].to_numpy()
        else:
            t = df.iloc[:, time_col].to_numpy()
            vin = df.iloc[:, input_col].to_numpy()
            vout = df.iloc[:, output_col].to_numpy()

        Ain, dAin, phi_in, dphi_in, Cin = fit_trace(t, vin, freq)
        Aout, dAout, phi_out, dphi_out, Cout = fit_trace(t, vout, freq)

        if Ain <= 0 or Aout <= 0:
            print(f"Skipping {file}: non-positive fitted amplitude")
            continue

        gain = Aout / Ain
        dgain = gain * np.sqrt((dAout / Aout)**2 + (dAin / Ain)**2)

        # Output phase relative to input, with automatic branch selection
        dphi = choose_phase_branch(
            phi_out, phi_in, freq,
            omega0_ref_guess, gamma_ref_guess
        )
        ddphi = np.sqrt(dphi_in**2 + dphi_out**2)

        frequencies.append(freq)
        gains.append(gain)
        gain_unc.append(dgain)
        phases.append(dphi)
        phase_unc.append(ddphi)

    except Exception as e:
        print(f"Skipping {file}: {e}")

frequencies = np.array(frequencies)
gains = np.array(gains)
gain_unc = np.array(gain_unc)
phases = np.array(phases)
phase_unc = np.array(phase_unc)

if len(frequencies) < 3:
    raise RuntimeError("Not enough valid datasets to perform Bode fit.")

# Sort by frequency
order = np.argsort(frequencies)
frequencies = frequencies[order]
gains = gains[order]
gain_unc = gain_unc[order]
phases = phases[order]
phase_unc = phase_unc[order]

# ============================================================
# FIT THE BODE MAGNITUDE CURVE
# ============================================================
f_peak_guess = frequencies[np.argmax(gains)]
omega0_guess = 2 * np.pi * f_peak_guess
gamma_guess = omega0_guess / 2
K_guess = np.max(gains) * gamma_guess * omega0_guess

p0 = [omega0_guess, gamma_guess, K_guess]
bounds = ([0, 0, 0], [np.inf, np.inf, np.inf])

popt_bode, pcov_bode = curve_fit(
    bode_amplitude,
    frequencies,
    gains,
    sigma=gain_unc,
    absolute_sigma=True,
    p0=p0,
    bounds=bounds,
    maxfev=20000
)

omega0_fit, gamma_fit, K_fit = popt_bode
domega0, dgamma, dK = np.sqrt(np.diag(pcov_bode))

# Derived resonance frequency in Hz
f0_fit = omega0_fit / (2 * np.pi)
df0_fit = domega0 / (2 * np.pi)

# Smooth fit curve
f_fit = np.logspace(np.log10(np.min(frequencies)), np.log10(np.max(frequencies)), 500)
gain_fit = bode_amplitude(f_fit, omega0_fit, gamma_fit, K_fit)

# ============================================================
# CONVERT TO dB
# ============================================================
gain_db = 20 * np.log10(gains)
gain_db_unc = 20 / np.log(10) * (gain_unc / gains)
gain_fit_db = 20 * np.log10(gain_fit)

# ============================================================
# PRINT RESULTS
# ============================================================
print("Bode fit results:")
print(f"omega0 = {omega0_fit:.3g} ± {domega0:.2g} rad s^-1")
print(f"f0     = {f0_fit:.3g} ± {df0_fit:.2g} Hz")
print(f"gamma  = {gamma_fit:.3g} ± {dgamma:.2g} s^-1")
print(f"K      = {K_fit:.3g} ± {dK:.2g}")

# ============================================================
# PLOTS
# ============================================================

# -------- Magnitude (linear) --------
plt.figure(figsize=(7, 5))
plt.errorbar(
    frequencies, gains, yerr=gain_unc,
    fmt='o', capsize=3, label='Data'
)
plt.semilogx(f_fit, gain_fit, label='Bode fit')
plt.xlabel("Frequency (Hz)")
plt.ylabel("Gain ($V_{out}/V_{in}$)")
plt.title("Bode Magnitude Plot")
plt.grid(True, which="both", ls="--", alpha=0.6)

textstr = (
    f"$f_0 = ({f0_fit:.3f} \\pm {df0_fit:.2g})\\,\\mathrm{{Hz}}$\n"
    f"$\\omega_0 = ({omega0_fit:.3f} \\pm {domega0:.2g})\\,\\mathrm{{rad\\,s^{{-1}}}}$\n"
    f"$\\gamma = ({gamma_fit:.3f} \\pm {dgamma:.2g})\\,\\mathrm{{s^{{-1}}}}$"
)
plt.text(
    0.03, 0.97, textstr,
    transform=plt.gca().transAxes,
    verticalalignment='top',
    bbox=dict(boxstyle='round', facecolor='white', alpha=0.85)
)

plt.legend()
plt.tight_layout()
plt.show()

# -------- Magnitude (dB) --------
plt.figure(figsize=(7, 5))
plt.errorbar(
    frequencies, gain_db, yerr=gain_db_unc,
    fmt='o', capsize=3, label='Data'
)
plt.semilogx(f_fit, gain_fit_db, label='Bode fit')
plt.xlabel("Frequency (Hz)")
plt.ylabel("Gain (dB)")
plt.title("Bode Magnitude Plot (dB)")
plt.grid(True, which="both", ls="--", alpha=0.6)

plt.text(
    0.03, 0.97, textstr,
    transform=plt.gca().transAxes,
    verticalalignment='top',
    bbox=dict(boxstyle='round', facecolor='white', alpha=0.85)
)

plt.legend()
plt.tight_layout()
plt.show()

# -------- Phase plot --------
phases_plot = np.unwrap(phases)
model_phase_plot = np.unwrap(bode_phase(f_fit, omega0_fit, gamma_fit))

plt.figure(figsize=(7, 5))
plt.errorbar(
    frequencies,
    np.degrees(phases_plot),
    yerr=np.degrees(phase_unc),
    fmt='o',
    capsize=3,
    label='Measured phase'
)
plt.semilogx(
    f_fit,
    np.degrees(model_phase_plot),
    label='Model phase'
)
plt.xlabel("Frequency (Hz)")
plt.ylabel("Phase (degrees)")
plt.title("Bode Phase Plot")
plt.grid(True, which="both", ls="--", alpha=0.6)
plt.legend()
plt.tight_layout()
plt.show()