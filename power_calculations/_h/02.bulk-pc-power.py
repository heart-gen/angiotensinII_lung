import numpy as np
from scipy.stats import norm
import matplotlib.pyplot as plt

def power_calc(n, mu=50, disp=0.3, fc=2.0, alpha=0.05, cv_between=0.4):
    """
    Estimate statistical power for detecting differential expression in bulk
    RNA-seq of isolated human lung pericytes. This uses NB model with
    inter-subject variability

    Parameters:
    - n (int): Number of biological replicates per group (COPD and control).
    - mu (float): Mean read count for the gene in the control group.
    - disp (float): Technical dispersion (from sequencing/technical noise).
    - fc (float): Fold change (treatment/control).
    - alpha (float): Significance threshold (e.g., 0.05).
    - cv_between (float): Coefficient of variation (biological/subject-level variability).

    Returns:
    - power (float): Estimated statistical power to detect the given fold change.
    """
    mu_treat = mu * fc
    # Biological variance
    # Moderate variability (e.g., homogeneous tissues): cv_between ≈ 0.3–0.4
    bio_var = (cv_between**2 * mu**2)
    bio_var_treat = (cv_between**2 * mu_treat**2)
    # Total variance includes technical (NB) and biological
    var_control = mu + disp * mu**2 + bio_var
    var_treat = mu_treat + disp * mu_treat**2 + bio_var_treat
    # Standard error of log2 fold change
    se_logfc = np.sqrt(var_control / (n * mu**2) + var_treat / (n * mu_treat**2))
    log_fc = np.log2(fc)
    # Power calculation
    z_alpha = norm.ppf(1 - alpha / 2)
    z_score = log_fc / se_logfc - z_alpha
    # Return power estimate
    return norm.cdf(z_score)


def plot_power_vs_dispersion(max_samples, dispersion_range, mu=50, fc=2.0,
                             cv_between=0.4, save_path='power_vs_dispersion.png'):
    """
    Plots and saves power vs. dispersion curves for different sample sizes in
    RNA-seq power analysis.

    Parameters:
        max_samples (int): Maximum sample size per group
        dispersion_range (np.ndarray): Array of dispersion values to test.
        mu (float): Mean expression.
        fc (float): Fold change.
        cv_between (float): Between-person coefficient of variation.
        save_path (str): File path to save the plot.
    """
    sample_sizes = range(2, max_samples+1)
    colors = plt.cm.viridis(np.linspace(0, 1, len(sample_sizes)))
    plt.figure(figsize=(8, 6))
    for n, color in zip(sample_sizes, colors):
        powers = [power_calc(n=n, disp=d) for d in dispersion_range]
        plt.plot(dispersion_range, powers, label=f'n = {n}', color=color)
    plt.axhline(0.8, color='gray', linestyle='--', label='80% power threshold')
    plt.xlabel('Dispersion (φ)', fontsize=18)
    plt.ylabel('Power', fontsize=18)
    plt.title('Power vs. Dispersion\nin Bulk RNA-seq (log2 FC=1)', fontsize=20)
    plt.xticks(fontsize=14); plt.yticks(fontsize=14)
    plt.legend(loc='lower left', title='Samples per Group')
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(save_path)
    plt.close()


def main():
    # Calculate required sample size
    max_samples = 10
    plot_power_vs_dispersion(max_samples,
                             dispersion_range=np.linspace(0.05, 0.5, 100))
