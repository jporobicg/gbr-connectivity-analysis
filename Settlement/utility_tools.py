#!/usr/bin/env python3
"""
Utility Functions for Larval Competency Modeling
================================================

This module contains all helper functions for data loading, preprocessing,
likelihood calculations, AICc computation, and statistical diagnostics.

All analysis logic should use functions from this module to ensure consistency.
"""

import pandas as pd
import numpy as np
import os
import re
from datetime import datetime
from typing import Tuple, Optional, Dict, Any, Callable
from scipy.integrate import quad
from scipy.optimize import minimize


# ============================================================================
# DATA LOADING AND PREPROCESSING
# ============================================================================

def get_family_name_mapping(abbreviations_path: Optional[str] = None) -> Dict[str, str]:
    """
    Get family name mapping to standardize names to match Species abbreviations.csv.
    
    This function creates a mapping from Settlement.csv family names to the
    standardized names used in Species abbreviations.csv.
    
    Parameters
    ----------
    abbreviations_path : str, optional
        Path to Species abbreviations.csv file. If None, uses default path.
        
    Returns
    -------
    Dict[str, str]
        Mapping from non-standard to standard family names
    """
    if abbreviations_path is None:
        # Default path relative to Settlement directory
        script_dir = os.path.dirname(os.path.abspath(__file__))
        abbreviations_path = os.path.join(
            script_dir, 'Randal_github', 'data', 'primary', 'Species abbreviations.csv'
        )
    
    # Standard mapping based on known inconsistencies
    # This maps Settlement.csv names to Species abbreviations.csv names
    mapping = {
        'Diploastreidae': 'Diploastraeidae',  # Missing 'a' in Settlement.csv
        'Lobophylliidae': 'Lobophyllidae',    # Extra 'i' in Settlement.csv
    }
    
    # If abbreviations file exists, verify and potentially extend mapping
    if os.path.exists(abbreviations_path):
        try:
            abbrev_df = pd.read_csv(abbreviations_path)
            standard_families = set(abbrev_df['Family'].dropna().unique())
            
            # Log any families in mapping that aren't in abbreviations file
            for old_name, new_name in mapping.items():
                if new_name not in standard_families:
                    print(f"WARNING: Mapped family '{new_name}' not found in Species abbreviations.csv")
        except Exception as e:
            print(f"WARNING: Could not read Species abbreviations file: {e}")
    
    return mapping


def standardize_family_names(
    df: pd.DataFrame,
    abbreviations_path: Optional[str] = None
) -> pd.DataFrame:
    """
    Standardize family names to match Species abbreviations.csv convention.
    
    This function applies the family name mapping to ensure consistency
    with the reference Species abbreviations.csv file.
    
    Parameters
    ----------
    df : pd.DataFrame
        DataFrame with 'Family' column to standardize
    abbreviations_path : str, optional
        Path to Species abbreviations.csv file
        
    Returns
    -------
    pd.DataFrame
        DataFrame with standardized family names
    """
    df = df.copy()
    
    if 'Family' not in df.columns:
        return df
    
    # Get mapping
    mapping = get_family_name_mapping(abbreviations_path)
    
    # Apply mapping
    df['Family'] = df['Family'].replace(mapping)
    
    # Log changes
    changes_made = []
    for old_name, new_name in mapping.items():
        if old_name in df['Family'].values:
            n_changed = (df['Family'] == old_name).sum()
            changes_made.append(f"{old_name} → {new_name} ({n_changed} rows)")
            df.loc[df['Family'] == old_name, 'Family'] = new_name
    
    if changes_made:
        print(f"Standardized family names:")
        for change in changes_made:
            print(f"  {change}")
    
    return df


def load_settlement_data(
    data_path: str,
    standardize_families: bool = True,
    abbreviations_path: Optional[str] = None,
    treatments: Optional[list] = None
) -> pd.DataFrame:
    """
    Load settlement data from CSV file.
    
    Parameters
    ----------
    data_path : str
        Path to the settlement data CSV file
    standardize_families : bool, optional
        If True, standardize family names to match Species abbreviations.csv, by default True
    abbreviations_path : str, optional
        Path to Species abbreviations.csv file. If None, uses default path.
    treatments : list, optional
        List of treatments to include. If None, includes all treatments. By default None.
        
    Returns
    -------
    pd.DataFrame
        DataFrame with columns: LarvalAge, NoSet, NoAlive, Family, etc.
        Filtered to remove invalid rows (missing Family, NoAlive <= 0)
        Family names standardized if standardize_families=True
        Filtered to specified treatments if treatments is provided
    """
    df = pd.read_csv(data_path)
    
    # Basic filtering
    df = df[df['Family'].notna()].copy()
    df = df[df['NoAlive'] > 0].copy()
    
    # Filter by treatment if specified
    if treatments is not None:
        if 'Treatment' not in df.columns:
            raise ValueError("Treatment column not found in data")
        df = df[df['Treatment'].isin(treatments)].copy()
        print(f"  Filtered to treatments: {treatments}")
        print(f"  Observations after treatment filter: {len(df)}")
    
    # Standardize family names to match Species abbreviations.csv
    if standardize_families:
        df = standardize_family_names(df, abbreviations_path)
    
    # Calculate observed proportion
    df['PropSettled'] = df['NoSet'] / df['NoAlive']
    
    return df


def load_monegetti_data(filepath: str) -> pd.DataFrame:
    """
    Load and process Monegetti's A. tenuis metamorphosis data.
    
    Parameters
    ----------
    filepath : str
        Path to Monegetti CSV file
        
    Returns
    -------
    pd.DataFrame
        Aggregated data by age with columns: LarvalAge, NoSet, NoAlive, PropSettled
    """
    df = pd.read_csv(filepath)
    
    # Rename columns for consistency
    df = df.rename(columns={
        'age (d)': 'LarvalAge',
        'meta': 'NoSet',
        'larvae': 'NoAlive'
    })
    
    # Calculate proportion metamorphosed (settled)
    df['PropSettled'] = df['NoSet'] / df['NoAlive']
    
    # Aggregate by age (sum across replicates)
    age_summary = df.groupby('LarvalAge').agg({
        'NoSet': 'sum',
        'NoAlive': 'sum'
    }).reset_index()
    age_summary['PropSettled'] = age_summary['NoSet'] / age_summary['NoAlive']
    age_summary = age_summary.sort_values('LarvalAge')
    
    return age_summary


# ============================================================================
# LIKELIHOOD CALCULATIONS
# ============================================================================

def neg_log_likelihood_binomial(
    params: np.ndarray,
    ages: np.ndarray,
    n_settled: np.ndarray,
    n_total: np.ndarray,
    model_func,
    phi: float = 1.0
) -> float:
    """
    Negative log-likelihood for binomial data with overdispersion.
    
    Parameters
    ----------
    params : np.ndarray
        Model parameters
    ages : np.ndarray
        Larval ages
    n_settled : np.ndarray
        Number settled (successes)
    n_total : np.ndarray
        Total number alive (trials)
    model_func : callable
        Competency function that takes (ages, *params) and returns predictions
    phi : float, optional
        Overdispersion parameter (1 = no overdispersion), by default 1.0
        
    Returns
    -------
    float
        Negative log-likelihood
    """
    # Predict competency
    pred_comp = model_func(ages, *params)
    pred_comp = np.clip(pred_comp, 1e-10, 1 - 1e-10)
    
    # Binomial log-likelihood
    ll = np.sum(
        n_settled * np.log(pred_comp) +
        (n_total - n_settled) * np.log(1 - pred_comp)
    )
    
    # Apply overdispersion correction
    ll = ll / phi
    
    return -ll


def neg_log_likelihood_monegetti(
    params: np.ndarray,
    ages: np.ndarray,
    n_settled: np.ndarray,
    n_total: np.ndarray,
    monegetti_func
) -> float:
    """
    Negative log-likelihood for Monegetti model with binomial data.
    
    Parameters
    ----------
    params : np.ndarray
        [a, b1, v1, b2, tc, Tcp] parameters
    ages : np.ndarray
        Larval ages
    n_settled : np.ndarray
        Number settled
    n_total : np.ndarray
        Total number alive
    monegetti_func : callable
        Monegetti competency function
        
    Returns
    -------
    float
        Negative log-likelihood
    """
    a, b1, v1, b2, tc, Tcp = params
    
    # Predict competency for each observation
    pred_comp = np.zeros_like(ages, dtype=float)
    for i, age in enumerate(ages):
        pred_comp[i] = monegetti_func(age, a, b1, v1, b2, tc, Tcp)
    
    # Clip to avoid log(0)
    pred_comp = np.clip(pred_comp, 1e-10, 1 - 1e-10)
    
    # Binomial log-likelihood
    ll = np.sum(
        n_settled * np.log(pred_comp) +
        (n_total - n_settled) * np.log(1 - pred_comp)
    )
    
    return -ll


# ============================================================================
# AICc COMPUTATION
# ============================================================================

def calculate_aicc(nll: float, n_params: int, n_obs: int) -> float:
    """
    Calculate AICc (Akaike Information Criterion corrected for small samples).
    
    Parameters
    ----------
    nll : float
        Negative log-likelihood
    n_params : int
        Number of parameters
    n_obs : int
        Number of observations
        
    Returns
    -------
    float
        AICc value (inf if n_obs <= n_params + 1)
    """
    aic = 2 * n_params + 2 * nll
    if n_obs > n_params + 1:
        correction = (2 * n_params * (n_params + 1)) / (n_obs - n_params - 1)
        return aic + correction
    else:
        return np.inf


# ============================================================================
# STATISTICAL DIAGNOSTICS
# ============================================================================

def estimate_overdispersion(
    df: pd.DataFrame,
    params: np.ndarray,
    model_func
) -> Tuple[float, np.ndarray]:
    """
    Estimate overdispersion parameter (phi) using Pearson residuals.
    
    phi > 1 indicates overdispersion
    
    Parameters
    ----------
    df : pd.DataFrame
        Data with columns: LarvalAge, NoSet, NoAlive
    params : np.ndarray
        Model parameters
    model_func : callable
        Competency function
        
    Returns
    -------
    Tuple[float, np.ndarray]
        (phi, pearson_residuals)
    """
    ages = df['LarvalAge'].values
    n_settled = df['NoSet'].values
    n_total = df['NoAlive'].values
    
    # Predicted proportions
    pred_prop = model_func(ages, *params)
    pred_prop = np.clip(pred_prop, 1e-10, 1 - 1e-10)
    
    # Observed proportions
    obs_prop = n_settled / n_total
    
    # Pearson residuals
    expected_var = pred_prop * (1 - pred_prop) / n_total
    pearson_resid = (obs_prop - pred_prop) / np.sqrt(expected_var)
    
    # Estimate phi (dispersion parameter)
    n_params = len(params)
    phi = np.sum(pearson_resid**2) / (len(df) - n_params)
    
    return phi, pearson_resid


def calculate_residuals(
    df: pd.DataFrame,
    params: np.ndarray,
    model_func
) -> Dict[str, np.ndarray]:
    """
    Calculate various residual types for model diagnostics.
    
    Parameters
    ----------
    df : pd.DataFrame
        Data with columns: LarvalAge, NoSet, NoAlive
    params : np.ndarray
        Model parameters
    model_func : callable
        Competency function
        
    Returns
    -------
    Dict[str, np.ndarray]
        Dictionary with 'pearson', 'deviance', 'raw' residuals
    """
    ages = df['LarvalAge'].values
    n_settled = df['NoSet'].values
    n_total = df['NoAlive'].values
    
    # Predicted proportions
    pred_prop = model_func(ages, *params)
    pred_prop = np.clip(pred_prop, 1e-10, 1 - 1e-10)
    
    # Observed proportions
    obs_prop = n_settled / n_total
    
    # Raw residuals
    raw_resid = obs_prop - pred_prop
    
    # Pearson residuals
    expected_var = pred_prop * (1 - pred_prop) / n_total
    pearson_resid = raw_resid / np.sqrt(expected_var)
    
    # Deviance residuals
    deviance_resid = np.sign(raw_resid) * np.sqrt(
        2 * (n_settled * np.log(n_settled / (n_total * pred_prop)) +
             (n_total - n_settled) * np.log((n_total - n_settled) /
                                           (n_total * (1 - pred_prop))))
    )
    deviance_resid = np.nan_to_num(deviance_resid, nan=0.0, posinf=0.0, neginf=0.0)
    
    return {
        'raw': raw_resid,
        'pearson': pearson_resid,
        'deviance': deviance_resid
    }


# ============================================================================
# FILE UTILITIES
# ============================================================================

def sanitize_filename(name: str) -> str:
    """
    Sanitize family name for use in filename.
    
    Removes or replaces special characters that are problematic in filenames.
    
    Parameters
    ----------
    name : str
        Family name to sanitize
        
    Returns
    -------
    str
        Sanitized name safe for use in filenames
    """
    # Replace spaces and special characters with underscores
    sanitized = re.sub(r'[^\w\s-]', '', str(name))
    sanitized = re.sub(r'[-\s]+', '_', sanitized)
    # Remove leading/trailing underscores
    sanitized = sanitized.strip('_')
    return sanitized


def ensure_directory(directory: str) -> None:
    """
    Ensure a directory exists, creating it if necessary.
    
    Parameters
    ----------
    directory : str
        Directory path
    """
    os.makedirs(directory, exist_ok=True)


# ============================================================================
# PARAMETER CONSTRAINTS
# ============================================================================

def get_monegetti_bounds(df: pd.DataFrame) -> list:
    """
    Get parameter bounds for Monegetti model based on data.
    
    Parameters
    ----------
    df : pd.DataFrame
        Data with columns: LarvalAge, NoSet
        
    Returns
    -------
    list
        List of (min, max) tuples for each parameter
    """
    ages = df['LarvalAge'].values
    
    # tc: must be before first settlement
    first_settlement = df[df['NoSet'] > 0]['LarvalAge'].min()
    tc_min = max(0, first_settlement - 2)
    tc_max = first_settlement
    
    # Tcp: must be after tc and before max age
    max_age = ages.max()
    Tcp_min = first_settlement
    Tcp_max = max_age
    
    bounds = [
        (0.1, 5.0),      # a
        (0.001, 1.0),    # b1
        (0.1, 5.0),      # v1
        (0.001, 1.0),    # b2
        (tc_min, tc_max),  # tc
        (Tcp_min, Tcp_max)  # Tcp
    ]
    
    return bounds


# ============================================================================
# TC50 CALCULATION
# ============================================================================

def calculate_tc50_monegetti(
    params: np.ndarray,
    max_age: float,
    monegetti_func,
    n_points: int = 100
) -> Tuple[Optional[float], float]:
    """
    Calculate TC50 (age at 50% of maximum competency) for Monegetti model.
    
    Parameters
    ----------
    params : np.ndarray
        [a, b1, v1, b2, tc, Tcp] parameters
    max_age : float
        Maximum age to test
    monegetti_func : callable
        Monegetti competency function
    n_points : int, optional
        Number of points to test, by default 100
        
    Returns
    -------
    Tuple[Optional[float], float]
        (TC50, max_competency). TC50 is None if max_competency < 0.1
    """
    a, b1, v1, b2, tc, Tcp = params
    
    # Test ages from tc to max_age
    test_ages = np.linspace(tc, max_age, n_points)
    pred_comp = np.array([
        monegetti_func(age, a, b1, v1, b2, tc, Tcp)
        for age in test_ages
    ])
    max_comp = pred_comp.max()
    
    if max_comp > 0.1:
        idx_50 = np.argmin(np.abs(pred_comp - 0.5 * max_comp))
        tc50 = test_ages[idx_50]
        return tc50, max_comp
    else:
        return None, max_comp


# ============================================================================
# COMPETENCY MODEL FUNCTIONS
# ============================================================================

def logistic_competency(age: np.ndarray, L: float, k: float, x0: float) -> np.ndarray:
    """
    Logistic (sigmoidal) competency model.
    
    Parameters
    ----------
    age : np.ndarray
        Larval ages in days
    L : float
        Maximum competency (asymptote)
    k : float
        Growth rate parameter
    x0 : float
        Inflection point (age at 50% of maximum)
        
    Returns
    -------
    np.ndarray
        Competency values
    """
    return L / (1 + np.exp(-k * (age - x0)))


def gompertz_competency(age: np.ndarray, a: float, b: float, c: float) -> np.ndarray:
    """
    Gompertz competency model.
    
    Parameters
    ----------
    age : np.ndarray
        Larval ages in days
    a : float
        Maximum competency (asymptote)
    b : float
        Displacement parameter
    c : float
        Growth rate parameter
        
    Returns
    -------
    np.ndarray
        Competency values
    """
    return a * np.exp(-b * np.exp(-c * age))


def weibull_competency(age: np.ndarray, tc: float, a: float, b: float, v: float) -> np.ndarray:
    """
    Weibull-based competency model.
    
    Parameters
    ----------
    age : np.ndarray
        Larval ages in days
    tc : float
        Precompetency period (onset age)
    a : float
        Maximum competency (asymptote)
    b : float
        Scale parameter
    v : float
        Shape parameter
        
    Returns
    -------
    np.ndarray
        Competency values
    """
    age = np.atleast_1d(age)
    competency = np.zeros_like(age, dtype=float)
    mask = age > tc
    if np.any(mask):
        competency[mask] = a * (1 - np.exp(-b * (age[mask] - tc)**v))
    return competency


def monegetti_competency(age: np.ndarray, a: float, b1: float, v1: float,
                         b2: float, tc: float, Tcp: float) -> np.ndarray:
    """
    Monegetti et al. piecewise Weibull-exponential competency model.
    
    Based on Moneghetti et al. (2019) "High-frequency sampling and piecewise
    models reshape dispersal kernels of a common reef coral"
    
    Model Description:
    - Precompetency period (t < tc): No settlement possible
    - Early competency (tc < t < Tcp): Weibull loss of competency
    - Late competency (t > Tcp): Exponential loss of competency
    
    Parameters
    ----------
    age : np.ndarray
        Larval age in days
    a : float
        Rate of acquisition of competency (when t > tc)
    b1 : float
        Weibull loss parameter (early period, scale)
    v1 : float
        Weibull loss parameter (early period, shape)
    b2 : float
        Exponential loss parameter (late period)
    tc : float
        Precompetency period (onset age)
    Tcp : float
        Change point between early and late phases
        
    Returns
    -------
    np.ndarray
        Probability of being competent at given age
    """
    age = np.atleast_1d(age)
    competency = np.zeros_like(age, dtype=float)
    
    for i, t in enumerate(age):
        if t < tc:
            # Precompetent
            competency[i] = 0.0
            
        elif t <= Tcp:
            # Early phase: integrate Weibull loss
            def integrand_early(tau):
                return a * np.exp(-a * (tau - tc)) * np.exp(-((b1 * (t - tau))**v1))
            
            try:
                result, _ = quad(integrand_early, tc, t, limit=50,
                               epsabs=1e-4, epsrel=1e-4)
                competency[i] = max(0, min(1, result))
            except Exception:
                competency[i] = 0.0
                
        else:  # t > Tcp
            # Late phase: integrate both phases
            # Part 1: Early phase up to Tcp
            def integrand_latter_1(tau):
                return (a * np.exp(-a * (tau - tc)) *
                       np.exp(-((b1 * (Tcp - tau))**v1)) *
                       np.exp(-b2 * (t - Tcp)))
            
            # Part 2: Late phase from Tcp to t
            def integrand_latter_2(tau):
                return a * np.exp(-a * (tau - tc)) * np.exp(-b2 * (t - tau))
            
            try:
                result1, _ = quad(integrand_latter_1, tc, Tcp, limit=50,
                                 epsabs=1e-4, epsrel=1e-4)
                result2, _ = quad(integrand_latter_2, Tcp, t, limit=50,
                                 epsabs=1e-4, epsrel=1e-4)
                competency[i] = max(0, min(1, result1 + result2))
            except Exception:
                competency[i] = 0.0
    
    return competency if len(competency) > 1 else competency[0]


# ============================================================================
# MODEL FITTING FUNCTIONS
# ============================================================================

def fit_logistic_binomial(
    df: pd.DataFrame,
    phi: float = 1.0
) -> Tuple[Optional[np.ndarray], float]:
    """
    Fit logistic model using binomial likelihood.
    
    Parameters
    ----------
    df : pd.DataFrame
        Data with columns: LarvalAge, NoSet, NoAlive
    phi : float, optional
        Overdispersion parameter, by default 1.0
        
    Returns
    -------
    Tuple[Optional[np.ndarray], float]
        (parameters, negative_log_likelihood) or (None, inf) if failed
    """
    ages = df['LarvalAge'].values
    n_settled = df['NoSet'].values
    n_total = df['NoAlive'].values
    
    def objective(params):
        L, k, x0 = params
        # Constraints
        if not (0.1 < L <= 1.0 and 0 < k < 10 and 0 < x0 < 50):
            return 1e10
        return neg_log_likelihood_binomial(
            params, ages, n_settled, n_total, logistic_competency, phi
        )
    
    # Try multiple starting points
    best_result = None
    best_nll = np.inf
    
    for L0 in [0.5, 0.7, 0.9]:
        for x0 in [3, 5, 7]:
            x0_init = [L0, 1.0, x0]
            try:
                result = minimize(
                    objective, x0_init, method='Nelder-Mead',
                    options={'maxiter': 5000}
                )
                if result.fun < best_nll and result.success:
                    best_nll = result.fun
                    best_result = result
            except Exception:
                continue
    
    if best_result is None or not best_result.success:
        return None, np.inf
    
    return best_result.x, best_result.fun


def fit_gompertz_binomial(
    df: pd.DataFrame,
    phi: float = 1.0
) -> Tuple[Optional[np.ndarray], float]:
    """
    Fit Gompertz model using binomial likelihood.
    
    Parameters
    ----------
    df : pd.DataFrame
        Data with columns: LarvalAge, NoSet, NoAlive
    phi : float, optional
        Overdispersion parameter, by default 1.0
        
    Returns
    -------
    Tuple[Optional[np.ndarray], float]
        (parameters, negative_log_likelihood) or (None, inf) if failed
    """
    ages = df['LarvalAge'].values
    n_settled = df['NoSet'].values
    n_total = df['NoAlive'].values
    
    def objective(params):
        a, b, c = params
        if not (0.1 < a <= 1.0 and 0 < b < 100 and 0 < c < 5):
            return 1e10
        return neg_log_likelihood_binomial(
            params, ages, n_settled, n_total, gompertz_competency, phi
        )
    
    best_result = None
    best_nll = np.inf
    
    for a0 in [0.5, 0.7, 0.9]:
        for c0 in [0.5, 1.0, 2.0]:
            x0_init = [a0, 1.0, c0]
            try:
                result = minimize(
                    objective, x0_init, method='Nelder-Mead',
                    options={'maxiter': 5000}
                )
                if result.fun < best_nll and result.success:
                    best_nll = result.fun
                    best_result = result
            except Exception:
                continue
    
    if best_result is None or not best_result.success:
        return None, np.inf
    
    return best_result.x, best_result.fun


def fit_weibull_binomial(
    df: pd.DataFrame,
    phi: float = 1.0
) -> Tuple[Optional[np.ndarray], float]:
    """
    Fit Weibull model using binomial likelihood.
    
    Parameters
    ----------
    df : pd.DataFrame
        Data with columns: LarvalAge, NoSet, NoAlive
    phi : float, optional
        Overdispersion parameter, by default 1.0
        
    Returns
    -------
    Tuple[Optional[np.ndarray], float]
        (parameters, negative_log_likelihood) or (None, inf) if failed
    """
    ages = df['LarvalAge'].values
    n_settled = df['NoSet'].values
    n_total = df['NoAlive'].values
    
    def objective(params):
        tc, a, b, v = params
        if not (0 < tc < 10 and 0.1 < a <= 2.0 and 0 < b < 5 and 0.1 < v < 5):
            return 1e10
        return neg_log_likelihood_binomial(
            params, ages, n_settled, n_total, weibull_competency, phi
        )
    
    best_result = None
    best_nll = np.inf
    
    for tc0 in [1, 2, 3]:
        for a0 in [0.5, 0.8]:
            x0_init = [tc0, a0, 0.1, 1.0]
            try:
                result = minimize(
                    objective, x0_init, method='Nelder-Mead',
                    options={'maxiter': 5000}
                )
                if result.fun < best_nll and result.success:
                    best_nll = result.fun
                    best_result = result
            except Exception:
                continue
    
    if best_result is None or not best_result.success:
        return None, np.inf
    
    return best_result.x, best_result.fun


def fit_monegetti_binomial(
    df: pd.DataFrame,
    monegetti_func: Callable
) -> Tuple[Optional[np.ndarray], float]:
    """
    Fit Monegetti piecewise model to settlement data.
    
    Parameters
    ----------
    df : pd.DataFrame
        Data with columns: LarvalAge, NoSet, NoAlive
    monegetti_func : Callable
        Monegetti competency function
        
    Returns
    -------
    Tuple[Optional[np.ndarray], float]
        (parameters, negative_log_likelihood) or (None, inf) if failed
    """
    ages = df['LarvalAge'].values
    n_settled = df['NoSet'].values
    n_total = df['NoAlive'].values
    
    # Get bounds
    bounds = get_monegetti_bounds(df)
    
    def objective(params):
        # Additional constraint: tc < Tcp
        if params[4] >= params[5]:
            return 1e10
        return neg_log_likelihood_monegetti(
            params, ages, n_settled, n_total, monegetti_func
        )
    
    # Try multiple starting points
    n_starts = 5
    best_result = None
    best_nll = np.inf
    
    # Create starting points
    tc_min, tc_max = bounds[4]
    Tcp_min, Tcp_max = bounds[5]
    
    starts = []
    for _ in range(n_starts):
        start = [
            np.random.uniform(0.5, 2.0),      # a
            np.random.uniform(0.01, 0.5),     # b1
            np.random.uniform(0.5, 2.0),      # v1
            np.random.uniform(0.01, 0.5),     # b2
            np.random.uniform(tc_min, tc_max),  # tc
            np.random.uniform(Tcp_min, Tcp_max)  # Tcp
        ]
        # Ensure tc < Tcp
        if start[4] >= start[5]:
            start[5] = start[4] + 1
        starts.append(start)
    
    for i, start in enumerate(starts):
        try:
            result = minimize(
                objective, start, method='L-BFGS-B',
                bounds=bounds,
                options={'maxiter': 50, 'ftol': 1e-3}
            )
            
            if result.fun < best_nll and result.success:
                best_nll = result.fun
                best_result = result
        except Exception:
            continue
    
    if best_result is not None and best_nll < 1e9:
        return best_result.x, best_nll
    else:
        return None, np.inf

