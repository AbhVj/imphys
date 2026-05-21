#!/usr/bin/env python3
"""
simple_fit_skeleton.py

Skeleton script for a simple 1D fit using iminuit.

"""
from scipy.stats import norm
import numpy as np
import matplotlib.pyplot as plt
from iminuit import Minuit

# ------------------------------
# 1. Made-up 1D data
# ------------------------------

# For simplicity, we'll create a 1D array of "measurements".
# In a real analysis, this might be something you have measured.
# Here we just generate some fixed numbers.

# Example: data drawn from a Gaussian with mean ~1.0 and sigma ~0.2 
rng = np.random.default_rng(42)
true_mean = 1.0
true_sigma = 0.2
n_data = 200
data = rng.normal(true_mean, true_sigma, size=n_data)

# 1D binning for the histogram
x_min, x_max = 0.0, 2.0
n_bins = 40
bin_edges = np.linspace(x_min, x_max, n_bins + 1)
bin_centres = 0.5 * (bin_edges[:-1] + bin_edges[1:])

# Build data histogram
data_hist, _ = np.histogram(data, bins=bin_edges)
print("[INFO] Data histogram counts:", data_hist)


# ------------------------------
# 2. Simple 1D model in parameter space
# ------------------------------

# We will use a simple Gaussian model with two parameters:
#   mu: mean
#   sigma: width
#
# The PDF (before normalisation) is:
#   f(x) = exp(-(x - mu)**2 / (2 * sigma**2))
#
# We will turn this into a binned expected count model for the histogram.


def unnormalised_gaussian(x, mu, sigma):
    return np.exp(-0.5 * ((x - mu) / sigma) ** 2)


# ------------------------------
# 3. Binned Poisson NLL for the 1D histogram
# ------------------------------

def nll_binned_factory(data_hist, bin_edges):
    """
    Create a binned Poisson NLL function for the histogram 'data_hist'.

    Parameters:
      data_hist : 1D array of observed counts in each bin
      bin_edges : edges of the bins (length = n_bins+1)

    The model will be a Gaussian in x with parameters mu and sigma,
    normalised so that the total expected counts equals the total
    number of observed counts.
    """
    data_counts = data_hist.astype(float)
    n_events = data_counts.sum()
    x_centres = 0.5 * (bin_edges[:-1] + bin_edges[1:])

    def nll(mu, sigma):
        # Protect against non-positive sigma
        if sigma <= 0:
            return 1e12

        # Unnormalised model at bin centres
        f_vals = unnormalised_gaussian(x_centres, mu, sigma)
        f_vals_sum = np.sum(f_vals)
        if f_vals_sum <= 0 or not np.isfinite(f_vals_sum):
            return 1e12

        # Normalise so that sum of expected counts = n_events
        pdf_vals = f_vals / f_vals_sum
        E = n_events * pdf_vals  # expected counts per bin
        O = data_counts

        # Poisson NLL = sum(E - O log E)
        mask = E > 0
        E_masked = E[mask]
        O_masked = O[mask]
        return np.sum(E_masked - O_masked * np.log(E_masked + 1e-12))

    return nll

nll = nll_binned_factory(data_hist, bin_edges)


# ------------------------------
# 4. Set up and run a simple Minuit fit
# ------------------------------

# Initial guesses for mu and sigma
initial_mu = 0.5
initial_sigma = 0.3

m = Minuit(
    nll,
    mu=initial_mu,
    sigma=initial_sigma,
)
m.errordef = Minuit.LIKELIHOOD  # for NLL
# [NEW] Apply lower bound on sigma
sigma_lower_bound = 0.17
m.limits["sigma"] = (sigma_lower_bound, None)

print("\n[INFO] Starting MIGRAD fit on simple 1D model...")
m.migrad()
print("[INFO] MIGRAD finished.")
print("[INFO] Running HESSE for errors...")
m.hesse()
print("[INFO] HESSE finished.\n")

print("Fit results (simple 1D Gaussian):")
print(f"  mu    = {m.values['mu']:.6f} ± {m.errors['mu']:.6f}")
print(f"  sigma = {m.values['sigma']:.6f} ± {m.errors['sigma']:.6f}")

# Covariance matrix for these two parameters
cov = m.covariance
cov_np = np.array(cov)
print("\nCovariance matrix (mu, sigma):")
print(cov_np)\
#5 Generate one toy histogram using the fitted parameters
def generate_toy_hist(mu, sigma, n_data, bin_edges, rng):
    """
    Generate one toy dataset of n_data points from a Gaussian with
    parameters (mu, sigma), histogram it using bin_edges, and return
    the toy histogram.

    Parameters
    ----------
    mu : float
        Gaussian mean.
    sigma : float
        Gaussian width.
    n_data : int
        Number of toy events to generate.
    bin_edges : np.ndarray
        Bin edges for the histogram.
    rng : np.random.Generator
        Random-number generator.

    Returns
    -------
    toy_hist : np.ndarray
        Histogram counts for the toy dataset.
    """
    if sigma <= 0:
        raise ValueError("sigma must be positive")

    # Generate one unbinned toy dataset
    toy_data = rng.normal(loc=mu, scale=sigma, size=n_data)

    # Histogram it with the same binning as the original data
    toy_hist, _ = np.histogram(toy_data, bins=bin_edges)

    return toy_hist
toy_hist = generate_toy_hist(
    mu=m.values["mu"],
    sigma=m.values["sigma"],
    n_data=n_data,
    bin_edges=bin_edges,
    rng=rng,
)

print("\n[INFO] Toy histogram counts:", toy_hist)
# ------------------------------
# 6. Run many toy experiments and refit each one
# ------------------------------

n_toys = 50

# Use the fitted values from the original data fit as the toy truth
fit_mu = m.values["mu"]
fit_sigma = m.values["sigma"]

# Arrays to store toy fit results
toy_mu = np.empty(n_toys)
toy_sigma = np.empty(n_toys)
toy_err_mu = np.empty(n_toys)
toy_err_sigma = np.empty(n_toys)

print(f"\n[INFO] Running {n_toys} toy experiments with sigma >= {sigma_lower_bound}...")

for i in range(n_toys):
    # Generate one toy histogram from the fitted model
    toy_hist = generate_toy_hist(
        mu=fit_mu,
        sigma=fit_sigma,
        n_data=n_data,
        bin_edges=bin_edges,
        rng=rng,
    )

    # Build a fresh NLL for this toy histogram
    toy_nll = nll_binned_factory(toy_hist, bin_edges)

    # Set up a Minuit fit for this toy
    toy_m = Minuit(
        toy_nll,
        mu=fit_mu,
        sigma=fit_sigma,
    )
    toy_m.errordef = Minuit.LIKELIHOOD
    # [NEW] Apply same sigma bound to toy fits
    toy_m.limits["sigma"] = (sigma_lower_bound, None)

    # Run the fit
    toy_m.migrad()
    toy_m.hesse()

    # Store fitted values and errors
    toy_mu[i] = toy_m.values["mu"]
    toy_sigma[i] = toy_m.values["sigma"]
    toy_err_mu[i] = toy_m.errors["mu"]
    toy_err_sigma[i] = toy_m.errors["sigma"]

print("[INFO] Toy study with bounded sigma finished.\n")

print("Stored toy fit results:")
print("toy_mu       =", toy_mu)
print("toy_sigma    =", toy_sigma)
print("toy_err_mu   =", toy_err_mu)
print("toy_err_sigma=", toy_err_sigma)
# ------------------------------
#7. Compute pull distributions
# ------------------------------

# Use the original fit results as the reference ("true") values
mu_fit = m.values["mu"]
sigma_fit = m.values["sigma"]

pull_mu = (toy_mu - mu_fit) / toy_err_mu
pull_sigma = (toy_sigma - sigma_fit) / toy_err_sigma

# Sample mean and sample standard deviation of the pull distributions
pull_mu_mean = np.mean(pull_mu)
pull_mu_std = np.std(pull_mu, ddof=1)

pull_sigma_mean = np.mean(pull_sigma)
pull_sigma_std = np.std(pull_sigma, ddof=1)

print("[INFO] Pull summary:")
print(f"  mu pull:    mean = {pull_mu_mean:.4f}, std = {pull_mu_std:.4f}")
print(f"  sigma pull: mean = {pull_sigma_mean:.4f}, std = {pull_sigma_std:.4f}")


# ------------------------------
# 8. Plot pull histograms
# ------------------------------

fig, axes = plt.subplots(1, 2, figsize=(12, 5))

# Pull plot for mu
axes[0].hist(pull_mu, bins=15, histtype="stepfilled", alpha=0.7, color="skyblue", edgecolor="black")
axes[0].axvline(pull_mu_mean, color="red", linestyle="--", linewidth=2, label=f"Mean = {pull_mu_mean:.3f}")
axes[0].set_xlabel("Pull of mu")
axes[0].set_ylabel("Number of toys")
axes[0].set_title("Pull distribution for mu")
axes[0].legend()

axes[0].text(
    0.05, 0.95,
    f"Mean = {pull_mu_mean:.3f}\nStd = {pull_mu_std:.3f}",
    transform=axes[0].transAxes,
    verticalalignment="top",
    bbox=dict(boxstyle="round", facecolor="white", alpha=0.8)
)

# Pull plot for sigma
axes[1].hist(pull_sigma, bins=15, histtype="stepfilled", alpha=0.7, color="lightgreen", edgecolor="black")
axes[1].axvline(pull_sigma_mean, color="red", linestyle="--", linewidth=2, label=f"Mean = {pull_sigma_mean:.3f}")
axes[1].set_xlabel("Pull of sigma")
axes[1].set_ylabel("Number of toys")
axes[1].set_title("Pull distribution for sigma")
axes[1].legend()

axes[1].text(
    0.05, 0.95,
    f"Mean = {pull_sigma_mean:.3f}\nStd = {pull_sigma_std:.3f}",
    transform=axes[1].transAxes,
    verticalalignment="top",
    bbox=dict(boxstyle="round", facecolor="white", alpha=0.8)
)

plt.tight_layout()
plt.savefig("pull_distributions.png", dpi=150)
plt.show()
# ------------------------------
# 9. Plot pull histogram for sigma
# ------------------------------

plt.figure(figsize=(7, 5))

plt.hist(
    pull_sigma,
    bins=15,
    histtype="stepfilled",
    alpha=0.7,
    color="salmon",
    edgecolor="black",
)

plt.axvline(
    pull_sigma_mean,
    color="red",
    linestyle="--",
    linewidth=2,
    label=f"Mean = {pull_sigma_mean:.3f}"
)

plt.xlabel("Pull of sigma")
plt.ylabel("Number of toys")
plt.title(f"Pull distribution for sigma with bound: sigma >= {sigma_lower_bound}")
plt.legend()

plt.text(
    0.05, 0.95,
    f"Mean = {pull_sigma_mean:.3f}\nStd = {pull_sigma_std:.3f}",
    transform=plt.gca().transAxes,
    verticalalignment="top",
    bbox=dict(boxstyle="round", facecolor="white", alpha=0.8)
)

plt.tight_layout()
plt.savefig("pull_sigma_bounded.png", dpi=150)
plt.show()
# ------------------------------
# 8. Fit Gaussians to the pull distributions
# ------------------------------

# Fit a Gaussian to each pull distribution
gauss_mu_mean, gauss_mu_sigma = norm.fit(pull_mu)
gauss_sigma_mean, gauss_sigma_sigma = norm.fit(pull_sigma)

# Approximate uncertainties on fitted Gaussian parameters
n_pull_mu = len(pull_mu)
n_pull_sigma = len(pull_sigma)

gauss_mu_mean_err = gauss_mu_sigma / np.sqrt(n_pull_mu)
gauss_mu_sigma_err = gauss_mu_sigma / np.sqrt(2 * n_pull_mu)

gauss_sigma_mean_err = gauss_sigma_sigma / np.sqrt(n_pull_sigma)
gauss_sigma_sigma_err = gauss_sigma_sigma / np.sqrt(2 * n_pull_sigma)

print("[INFO] Gaussian fit to pull distributions:")
print(
    f"  mu pull Gaussian: mean = {gauss_mu_mean:.4f} ± {gauss_mu_mean_err:.4f}, "
    f"sigma = {gauss_mu_sigma:.4f} ± {gauss_mu_sigma_err:.4f}"
)
print(
    f"  sigma pull Gaussian: mean = {gauss_sigma_mean:.4f} ± {gauss_sigma_mean_err:.4f}, "
    f"sigma = {gauss_sigma_sigma:.4f} ± {gauss_sigma_sigma_err:.4f}"
)
# ------------------------------
# 9. Plot pull histograms with Gaussian overlays
# ------------------------------

fig, axes = plt.subplots(1, 2, figsize=(12, 5))

# Common x-ranges for smooth Gaussian curves
x_mu = np.linspace(np.min(pull_mu), np.max(pull_mu), 400)
x_sigma = np.linspace(np.min(pull_sigma), np.max(pull_sigma), 400)

# ---- mu pull ----
counts_mu, bins_mu, _ = axes[0].hist(
    pull_mu,
    bins=15,
    histtype="stepfilled",
    alpha=0.7,
    color="skyblue",
    edgecolor="black",
    density=True,
    label="Toy pulls"
)

axes[0].plot(
    x_mu,
    norm.pdf(x_mu, gauss_mu_mean, gauss_mu_sigma),
    color="darkblue",
    linewidth=2,
    label="Gaussian fit"
)

axes[0].axvline(0.0, color="black", linestyle=":", linewidth=1.5, label="Zero pull")
axes[0].axvline(gauss_mu_mean, color="red", linestyle="--", linewidth=2, label=f"Fit mean = {gauss_mu_mean:.3f}")

axes[0].set_xlabel("Pull of mu")
axes[0].set_ylabel("Density")
axes[0].set_title("Pull distribution for mu")

axes[0].text(
    0.05, 0.95,
    (
        f"Gaussian fit:\n"
        f"mean = {gauss_mu_mean:.3f} ± {gauss_mu_mean_err:.3f}\n"
        f"sigma = {gauss_mu_sigma:.3f} ± {gauss_mu_sigma_err:.3f}"
    ),
    transform=axes[0].transAxes,
    verticalalignment="top",
    bbox=dict(boxstyle="round", facecolor="white", alpha=0.8)
)

axes[0].legend()

# ---- sigma pull ----
counts_sigma, bins_sigma, _ = axes[1].hist(
    pull_sigma,
    bins=15,
    histtype="stepfilled",
    alpha=0.7,
    color="salmon",
    edgecolor="black",
    density=True,
    label="Toy pulls"
)

axes[1].plot(
    x_sigma,
    norm.pdf(x_sigma, gauss_sigma_mean, gauss_sigma_sigma),
    color="darkred",
    linewidth=2,
    label="Gaussian fit"
)

axes[1].axvline(0.0, color="black", linestyle=":", linewidth=1.5, label="Zero pull")
axes[1].axvline(gauss_sigma_mean, color="red", linestyle="--", linewidth=2, label=f"Fit mean = {gauss_sigma_mean:.3f}")

axes[1].set_xlabel("Pull of sigma")
axes[1].set_ylabel("Density")
axes[1].set_title("Pull distribution for sigma")

axes[1].text(
    0.05, 0.95,
    (
        f"Gaussian fit:\n"
        f"mean = {gauss_sigma_mean:.3f} ± {gauss_sigma_mean_err:.3f}\n"
        f"sigma = {gauss_sigma_sigma:.3f} ± {gauss_sigma_sigma_err:.3f}"
    ),
    transform=axes[1].transAxes,
    verticalalignment="top",
    bbox=dict(boxstyle="round", facecolor="white", alpha=0.8)
)

axes[1].legend()

plt.tight_layout()
plt.savefig("pull_distributions_with_gaussian_fit.png", dpi=150)
plt.show()