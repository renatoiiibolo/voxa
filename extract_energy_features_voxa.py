"""
extract_energy_features_voxa.py — Local energy feature extraction for the VA extension.

Reads TOPAS-nBio DSB coordinate files and extracts per-DSB energy features
used to calibrate and apply the Voxel-Aware (VA) extension in steps 9–11.

The primary feature is E_local: energy deposited in the DSB's own chromatin
voxel. Additional neighborhood and gradient features are also computed but
E_local (via its z-score) is what delta_f acts on in the VA formulation:

    f_direct(i) = f_dir(x) + delta_f * E_zscore(i)

Usage:
    # Calibration arm (~2500 DSBs per particle, used in step 9)
    python extract_energy_features_voxa.py --mode calibration \
        --data-dir ./calibration_data \
        --output-dir ./voxa_features_output_calibration

    # Validation arm (~400 DSBs per particle, used in steps 11–12)
    python extract_energy_features_voxa.py --mode validation \
        --data-dir ./validation_data \
        --output-dir ./voxa_features_output_validation

Outputs per run:
    {particle}_{mode}_dsb_energy_features.csv
    {particle}_{mode}_energy_statistics.json
    all_particles_{mode}_energy_features.csv
    combined_{mode}_statistics.json
"""

================================================================================
ENERGY FEATURE EXTRACTION FOR VOxA MODEL
================================================================================

Version:  VOxA v1.0 (Variable Oxygen-dependent Amorphous Track Model)
Date:     February 2026

DESCRIPTION
-----------
This script extracts local energy deposition features from TOPAS-nBio Monte
Carlo simulations for use with the Voxel-Aware (VA) extension of the VOxA model.

The VOxA model has two parts:
  1. OM (Oxygen Model): Population-level OER predictions (Steps 1-8)
  2. VA (Voxel-Aware): DSB-level P_DSB predictions using local energy (Steps 9-11)

This script supports the VA workflow:
  - Step 9:  Calibration using 2500-DSB datasets
  - Step 11: Validation using 400-DSB datasets

KEY FEATURES EXTRACTED
----------------------
- E_local: Energy in the DSB's voxel (PRIMARY feature for VA model)
- E_total_3x3x3: Total energy in 3×3×3 neighborhood
- E_mean_3x3x3: Mean energy in neighborhood
- Gradient features: Local energy gradient magnitude
- Core fraction: E_local / E_total (track core indicator)

VA MODEL FORMULATION
--------------------
The VA extension modulates the direct damage fraction for each DSB:

  f_direct(i) = p1_base + δf × E_zscore(i)

where:
  - p1_base: Particle-specific baseline direct fraction (from OM calibration)
  - δf: Energy sensitivity parameter (calibrated per particle in Step 9)
  - E_zscore(i): Standardized local energy for DSB i

OUTPUT FILES
------------
- {particle}_{mode}_dsb_energy_features.csv: Per-DSB features
- {particle}_{mode}_energy_statistics.json: Summary statistics
- all_particles_{mode}_energy_features.csv: Combined file for R calibration
- combined_{mode}_statistics.json: Combined statistics

USAGE
-----
  # Extract features from 2500-DSB calibration datasets (Step 9)
  python extract_energy_features_voxa.py --mode calibration \\
      --data-dir ./calibration_data \\
      --output-dir ./voxa_features_output_calibration

  # Extract features from 400-DSB validation datasets (Step 11)
  python extract_energy_features_voxa.py --mode validation \\
      --data-dir ./validation_data \\
      --output-dir ./voxa_features_output_validation

================================================================================
"""

import argparse
import numpy as np
import pandas as pd
from pathlib import Path
from dataclasses import dataclass, asdict
from typing import Tuple, Dict, List
import json
import logging
from datetime import datetime

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


# VOxA MODEL REFERENCE PARAMETERS

VOXA_REFERENCE = {
    "model_version": "VOxA v1.0 (Variable Oxygen-dependent Amorphous Track Model)",
    "model_date": "2026-02",
    
    # Oxygen kinetics from OM calibration (Step 3)
    "K_fix": 0.1593,      # % O2
    "K_repair": 0.2119,   # % O2
    "OER_max_retention": 3.32,
    "OER_max_survival": 2.93,
    "conversion_factor": 1.20,  # retention to survival
    
    # Case fractions from DSB combinatorics (Sakata 2019)
    "case_fractions": {
        "p1_low": 0.04,   # d² = purely direct DSBs
        "p2_low": 0.32,   # 2di = hybrid DSBs
        "p3_low": 0.64,   # i² = purely indirect DSBs
        "p1_high": 0.68,  # direct at high LET (Hirayama 2009)
    },
    
    # Physical bounds for f_direct
    "f_min": 0.02,  # minimum direct fraction
    "f_max": 0.68,  # maximum direct fraction (Hirayama 2009)
    
    # Calibration oxygen level for VA
    "calibration_O2": 0.21,  # % O2 (moderate hypoxia)
    
    # Particle-specific parameters (from VA calibration, Step 9)
    # These will be updated after running Step 9
    "particles": {
        "electron": {
            "delta_f": 0.00357,
            "delta_f_ci_low": 0.00348,
            "delta_f_ci_high": 0.00486,
            "P_DSB_cv": 0.36,
            "n_calibration": 2515,
            "p1_base": 0.0401,
            "p2_base": 0.3199,
            "p3_base": 0.6400,
            "LET_nominal": 0.2,  # keV/μm
        },
        "proton": {
            "delta_f": 0.02588,
            "delta_f_ci_low": 0.02431,
            "delta_f_ci_high": 0.02745,
            "P_DSB_cv": 2.58,
            "n_calibration": 2491,
            "p1_base": 0.0953,
            "p2_base": 0.2649,
            "p3_base": 0.6398,
            "LET_nominal": 4.6,
        },
        "carbon": {
            "delta_f": 0.08753,
            "delta_f_ci_low": 0.08155,
            "delta_f_ci_high": 0.09365,
            "P_DSB_cv": 8.08,
            "n_calibration": 2519,
            "p1_base": 0.3879,
            "p2_base": 0.0000,
            "p3_base": 0.6121,
            "LET_nominal": 40.9,
        },
    }
}


# GEOMETRY CONSTANTS

@dataclass
class SimulationGeometry:
    """Complete geometry parameters for TOPAS-nBio simulation."""
    nucleus_radius_um: float = 4.65
    box_half_length_um: float = 9.28055
    n_energy_voxels: int = 60
    n_dose_voxels: int = 30

    @property
    def box_full_length_um(self) -> float:
        return 2 * self.box_half_length_um

    @property
    def energy_voxel_size_um(self) -> float:
        return self.box_full_length_um / self.n_energy_voxels

    @property
    def dose_voxel_size_um(self) -> float:
        return self.box_full_length_um / self.n_dose_voxels

    def validate_against_topas_header(self, bin_size_cm: float, grid_type: str = "energy") -> bool:
        """Validate geometry against TOPAS file header."""
        bin_size_um = bin_size_cm * 1e4
        expected = self.energy_voxel_size_um if grid_type == "energy" else self.dose_voxel_size_um
        match = np.isclose(bin_size_um, expected, rtol=0.01)
        status = "✓ PASSED" if match else "⚠ WARNING"
        logger.info(f"Geometry validation {status} ({grid_type}):")
        logger.info(f"  TOPAS bin size: {bin_size_um:.6f} µm")
        logger.info(f"  Expected:       {expected:.6f} µm")
        return match

    def print_summary(self) -> None:
        """Print geometry configuration summary."""
        print("\n" + "=" * 70)
        print("SIMULATION GEOMETRY")
        print("=" * 70)
        print(f"\nParallelWorldBox:")
        print(f"  Half-length: {self.box_half_length_um} µm")
        print(f"  Full length: {self.box_full_length_um:.4f} µm")
        print(f"  Bounds: [{-self.box_half_length_um}, +{self.box_half_length_um}] µm")
        print(f"\nNucleus (centered in box):")
        print(f"  Radius: {self.nucleus_radius_um} µm")
        print(f"\nEnergy Grid ({self.n_energy_voxels}³):")
        print(f"  Voxel size: {self.energy_voxel_size_um:.6f} µm")
        nucleus_voxels = int(2 * self.nucleus_radius_um / self.energy_voxel_size_um)
        print(f"  ~{nucleus_voxels} voxels across nucleus diameter")
        print("=" * 70)

    def to_dict(self) -> Dict:
        """Export geometry as dictionary."""
        return {
            'box_half_length_um': self.box_half_length_um,
            'box_full_length_um': self.box_full_length_um,
            'nucleus_radius_um': self.nucleus_radius_um,
            'n_energy_voxels': self.n_energy_voxels,
            'energy_voxel_size_um': self.energy_voxel_size_um,
        }


# COORDINATE TO VOXEL MAPPING

def coord_to_energy_voxel(
    x: float, y: float, z: float,
    geom: SimulationGeometry = None
) -> Tuple[int, int, int]:
    """Map 3D coordinate (µm) to energy voxel indices."""
    if geom is None:
        geom = SimulationGeometry()

    H = geom.box_half_length_um
    n = geom.n_energy_voxels
    delta = geom.energy_voxel_size_um

    # Clip to box bounds
    x = np.clip(x, -H + 1e-9, H - 1e-9)
    y = np.clip(y, -H + 1e-9, H - 1e-9)
    z = np.clip(z, -H + 1e-9, H - 1e-9)

    # Convert to voxel indices
    i = int(np.floor((x + H) / delta))
    j = int(np.floor((y + H) / delta))
    k = int(np.floor((z + H) / delta))

    return (np.clip(i, 0, n-1), np.clip(j, 0, n-1), np.clip(k, 0, n-1))


def get_neighboring_voxels(
    i_center: int, j_center: int, k_center: int,
    radius: int = 1,
    n_voxels: int = 60
) -> List[Tuple[int, int, int]]:
    """Get indices of neighboring voxels within a given radius."""
    neighbors = []
    for di in range(-radius, radius + 1):
        for dj in range(-radius, radius + 1):
            for dk in range(-radius, radius + 1):
                ni, nj, nk = i_center + di, j_center + dj, k_center + dk
                if 0 <= ni < n_voxels and 0 <= nj < n_voxels and 0 <= nk < n_voxels:
                    neighbors.append((ni, nj, nk))
    return neighbors


# FILE PARSING

def detect_delimiter(filepath: Path) -> str:
    """Detect file delimiter from first few lines."""
    with open(filepath, 'r') as f:
        lines = [f.readline() for _ in range(5)]
    
    # Skip comment lines
    data_lines = [l for l in lines if l.strip() and not l.strip().startswith('#')]
    
    if not data_lines:
        return ','
    
    first_data = data_lines[0]
    
    if '\t' in first_data and first_data.count('\t') >= 2:
        return '\t'
    elif ',' in first_data:
        return ','
    else:
        return r'\s+'


def parse_dsb_file(
    filepath: Path,
    coord_columns: Tuple[str, str, str] = ('x_um', 'y_um', 'z_um')
) -> pd.DataFrame:
    """Parse DSB coordinates CSV file."""
    filepath = Path(filepath)
    logger.info(f"Parsing DSB file: {filepath.name}")

    sep = detect_delimiter(filepath)
    sep_name = "TAB" if sep == '\t' else ("COMMA" if sep == ',' else "WHITESPACE")
    logger.info(f"  Detected delimiter: {sep_name}")

    df = pd.read_csv(filepath, sep=sep)

    logger.info(f"  Columns found: {list(df.columns)[:6]}{'...' if len(df.columns) > 6 else ''}")
    logger.info(f"  Total DSBs: {len(df)}")

    # Verify coordinate columns
    x_col, y_col, z_col = coord_columns
    missing = [c for c in coord_columns if c not in df.columns]
    if missing:
        # Try common alternatives
        alt_mappings = {
            'x_um': ['x', 'X', 'x_position', 'pos_x'],
            'y_um': ['y', 'Y', 'y_position', 'pos_y'],
            'z_um': ['z', 'Z', 'z_position', 'pos_z'],
        }
        for col in missing:
            for alt in alt_mappings.get(col, []):
                if alt in df.columns:
                    df = df.rename(columns={alt: col})
                    logger.info(f"  Mapped column '{alt}' -> '{col}'")
                    break
        
        # Check again
        missing = [c for c in coord_columns if c not in df.columns]
        if missing:
            raise ValueError(f"Missing coordinate columns: {missing}\nAvailable: {list(df.columns)}")

    logger.info(f"  X range: [{df[x_col].min():.3f}, {df[x_col].max():.3f}] µm")
    logger.info(f"  Y range: [{df[y_col].min():.3f}, {df[y_col].max():.3f}] µm")
    logger.info(f"  Z range: [{df[z_col].min():.3f}, {df[z_col].max():.3f}] µm")

    return df


def parse_topas_energy_file(
    filepath: Path,
    expected_shape: Tuple[int, int, int] = (60, 60, 60)
) -> Tuple[np.ndarray, Dict]:
    """
    Parse TOPAS-nBio EnergyDeposit file.
    
    Handles multiple formats:
    - Tab-separated: i<tab>j<tab>k<tab>value
    - Comma-separated: i, j, k, value
    - Space-separated: i j k value
    """
    filepath = Path(filepath)
    logger.info(f"Parsing TOPAS energy file: {filepath.name}")

    energy_grid = np.zeros(expected_shape)
    metadata = {
        'filepath': str(filepath),
        'header_lines': [],
        'bin_size_cm': None,
        'n_bins': expected_shape,
    }

    n_data_lines = 0
    n_nonzero = 0
    parse_errors = 0

    # Detect delimiter
    detected_delimiter = None
    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('#') or not line:
                continue
            if '\t' in line:
                detected_delimiter = 'tab'
            elif ', ' in line:
                detected_delimiter = 'comma-space'
            elif ',' in line:
                detected_delimiter = 'comma'
            else:
                detected_delimiter = 'space'
            break

    logger.info(f"  Detected data format: {detected_delimiter}")

    with open(filepath, 'r') as f:
        for line_num, line in enumerate(f):
            line_stripped = line.strip()

            if not line_stripped:
                continue

            if line_stripped.startswith('#'):
                metadata['header_lines'].append(line_stripped)
                # Extract bin size from header
                if 'bins of' in line_stripped.lower():
                    try:
                        parts = line_stripped.split('bins of')
                        if len(parts) == 2:
                            size_str = parts[1].strip().split()[0]
                            metadata['bin_size_cm'] = float(size_str)
                    except (ValueError, IndexError):
                        pass
                continue

            # Parse data line
            try:
                if detected_delimiter in ['comma-space', 'comma']:
                    parts = [p.strip() for p in line_stripped.split(',')]
                elif detected_delimiter == 'tab':
                    parts = line_stripped.split('\t')
                else:
                    parts = line_stripped.split()

                if len(parts) >= 4:
                    i = int(parts[0])
                    j = int(parts[1])
                    k = int(parts[2])
                    value = float(parts[3])

                    if (0 <= i < expected_shape[0] and
                        0 <= j < expected_shape[1] and
                        0 <= k < expected_shape[2]):
                        energy_grid[i, j, k] = value
                        n_data_lines += 1
                        if value > 0:
                            n_nonzero += 1
                else:
                    parse_errors += 1

            except (ValueError, IndexError) as e:
                parse_errors += 1
                if parse_errors <= 3:
                    logger.warning(f"  Line {line_num}: parse error: {e}")
                continue

    n_total = np.prod(expected_shape)
    total_energy = np.sum(energy_grid)
    max_energy = np.max(energy_grid)

    metadata['n_data_lines'] = n_data_lines
    metadata['n_nonzero_voxels'] = n_nonzero
    metadata['total_energy_MeV'] = float(total_energy)
    metadata['max_voxel_energy_MeV'] = float(max_energy)
    metadata['fraction_nonzero'] = n_nonzero / n_total if n_total > 0 else 0
    metadata['parse_errors'] = parse_errors

    logger.info(f"  Data lines parsed: {n_data_lines:,}")
    logger.info(f"  Non-zero voxels: {n_nonzero:,} ({100*metadata['fraction_nonzero']:.2f}%)")
    logger.info(f"  Total energy: {total_energy:.6f} MeV")
    logger.info(f"  Max voxel energy: {max_energy:.6e} MeV")

    if metadata['bin_size_cm']:
        bin_um = metadata['bin_size_cm'] * 1e4
        logger.info(f"  Bin size from header: {bin_um:.4f} µm")

    return energy_grid, metadata


# ENERGY FEATURE COMPUTATION

@dataclass
class DSBEnergyFeatures:
    """Energy features for a single DSB location."""
    dsb_id: int
    x: float
    y: float
    z: float
    complexity: str = ""
    i_energy: int = 0
    j_energy: int = 0
    k_energy: int = 0
    
    # Primary feature for VA model
    E_local: float = 0.0
    
    # Neighborhood features
    E_total_3x3x3: float = 0.0
    E_mean_3x3x3: float = 0.0
    E_max_3x3x3: float = 0.0
    E_std_3x3x3: float = 0.0
    n_active_3x3x3: int = 0
    
    # Gradient features
    grad_x: float = 0.0
    grad_y: float = 0.0
    grad_z: float = 0.0
    grad_mag: float = 0.0
    
    # Derived features
    core_frac: float = 0.0
    is_local_max: bool = False
    
    # Metadata
    run_id: str = ""
    particle: str = ""
    dataset_type: str = ""  # "calibration" or "validation"


def compute_gradient(
    energy_grid: np.ndarray,
    center: Tuple[int, int, int]
) -> Tuple[float, float, float]:
    """Compute energy gradient at center voxel using central differences."""
    i, j, k = center
    n = energy_grid.shape[0]

    grad_x = grad_y = grad_z = 0.0

    if 0 < i < n - 1:
        grad_x = (energy_grid[i+1, j, k] - energy_grid[i-1, j, k]) / 2
    elif i == 0 and n > 1:
        grad_x = energy_grid[i+1, j, k] - energy_grid[i, j, k]
    elif i == n - 1 and n > 1:
        grad_x = energy_grid[i, j, k] - energy_grid[i-1, j, k]

    if 0 < j < n - 1:
        grad_y = (energy_grid[i, j+1, k] - energy_grid[i, j-1, k]) / 2
    elif j == 0 and n > 1:
        grad_y = energy_grid[i, j+1, k] - energy_grid[i, j, k]
    elif j == n - 1 and n > 1:
        grad_y = energy_grid[i, j, k] - energy_grid[i, j-1, k]

    if 0 < k < n - 1:
        grad_z = (energy_grid[i, j, k+1] - energy_grid[i, j, k-1]) / 2
    elif k == 0 and n > 1:
        grad_z = energy_grid[i, j, k+1] - energy_grid[i, j, k]
    elif k == n - 1 and n > 1:
        grad_z = energy_grid[i, j, k] - energy_grid[i, j, k-1]

    return float(grad_x), float(grad_y), float(grad_z)


def extract_dsb_energy_features(
    dsb_df: pd.DataFrame,
    energy_grid: np.ndarray,
    coord_columns: Tuple[str, str, str] = ('x_um', 'y_um', 'z_um'),
    complexity_column: str = 'complexity',
    geom: SimulationGeometry = None,
    neighborhood_radius: int = 1,
    particle: str = "",
    run_id: str = "",
    dataset_type: str = "calibration"
) -> List[DSBEnergyFeatures]:
    """Extract energy features for all DSB locations."""
    if geom is None:
        geom = SimulationGeometry()

    x_col, y_col, z_col = coord_columns
    features_list = []
    n_voxels = geom.n_energy_voxels

    for idx, row in dsb_df.iterrows():
        x, y, z = row[x_col], row[y_col], row[z_col]

        complexity = ""
        if complexity_column in dsb_df.columns:
            complexity = str(row[complexity_column])

        # Map to voxel
        i_e, j_e, k_e = coord_to_energy_voxel(x, y, z, geom)
        E_local = energy_grid[i_e, j_e, k_e]

        # Neighborhood analysis
        neighbors = get_neighboring_voxels(i_e, j_e, k_e, neighborhood_radius, n_voxels)
        neighbor_energies = np.array([energy_grid[ni, nj, nk] for ni, nj, nk in neighbors])

        E_total_3x3x3 = float(np.sum(neighbor_energies))
        E_mean_3x3x3 = float(np.mean(neighbor_energies))
        E_max_3x3x3 = float(np.max(neighbor_energies))
        E_std_3x3x3 = float(np.std(neighbor_energies))
        n_active_3x3x3 = int(np.sum(neighbor_energies > 0))

        # Gradient
        grad_x, grad_y, grad_z = compute_gradient(energy_grid, (i_e, j_e, k_e))
        grad_mag = float(np.sqrt(grad_x**2 + grad_y**2 + grad_z**2))

        # Derived features
        eps = 1e-12
        core_frac = E_local / (E_total_3x3x3 + eps) if E_total_3x3x3 > 0 else 0.0
        is_local_max = (E_local >= E_max_3x3x3) and (E_local > 0)

        features = DSBEnergyFeatures(
            dsb_id=int(idx),
            x=float(x), y=float(y), z=float(z),
            complexity=complexity,
            i_energy=i_e, j_energy=j_e, k_energy=k_e,
            E_local=E_local,
            E_total_3x3x3=E_total_3x3x3,
            E_mean_3x3x3=E_mean_3x3x3,
            E_max_3x3x3=E_max_3x3x3,
            E_std_3x3x3=E_std_3x3x3,
            n_active_3x3x3=n_active_3x3x3,
            grad_x=grad_x, grad_y=grad_y, grad_z=grad_z,
            grad_mag=grad_mag,
            core_frac=float(core_frac),
            is_local_max=is_local_max,
            run_id=run_id,
            particle=particle,
            dataset_type=dataset_type
        )
        features_list.append(features)

    return features_list


def features_to_dataframe(features: List[DSBEnergyFeatures]) -> pd.DataFrame:
    """Convert features list to DataFrame."""
    return pd.DataFrame([asdict(f) for f in features])


# VOxA VA MODEL PREDICTIONS

def compute_va_predictions(
    df: pd.DataFrame,
    particle: str,
    O2_levels: List[float] = None,
    delta_f: float = None
) -> pd.DataFrame:
    """
    Compute VOxA Voxel-Aware predictions for extracted features.
    
    This applies the calibrated δf values to predict P_DSB for each DSB.
    
    Parameters:
    -----------
    df : pd.DataFrame
        DataFrame with E_local column
    particle : str
        Particle type (electron, proton, carbon)
    O2_levels : List[float]
        Oxygen levels to compute predictions for (% O2)
    delta_f : float, optional
        Override δf value (for sensitivity analysis)
    
    Returns:
    --------
    pd.DataFrame
        Input DataFrame with added P_DSB columns
    """
    if O2_levels is None:
        O2_levels = [21.0, 2.1, 0.21, 0.021, 0.001]
    
    particle_lower = particle.lower()
    if particle_lower not in VOXA_REFERENCE['particles']:
        logger.warning(f"Unknown particle '{particle}', using electron parameters")
        particle_lower = 'electron'
    
    params = VOXA_REFERENCE['particles'][particle_lower]
    K_fix = VOXA_REFERENCE['K_fix']
    K_repair = VOXA_REFERENCE['K_repair']
    f_min = VOXA_REFERENCE['f_min']
    f_max = VOXA_REFERENCE['f_max']
    
    # Use provided delta_f or calibrated value
    if delta_f is None:
        delta_f = params['delta_f']
    
    # Get E_local and compute z-scores
    E_local = df['E_local'].values
    E_mean = np.mean(E_local)
    E_std = np.std(E_local)
    
    if E_std < 1e-10:
        E_zscore = np.zeros_like(E_local)
        logger.warning("E_local has zero variance - all z-scores set to 0")
    else:
        E_zscore = (E_local - E_mean) / E_std
    
    # Base case fractions from calibration
    p1_base = params['p1_base']
    p2_base = params['p2_base']
    p3_base = params['p3_base']
    
    # Compute f_direct for each DSB (Equation 7 from formulation)
    f_direct = np.clip(p1_base + delta_f * E_zscore, f_min, f_max)
    
    # Store intermediate values
    df = df.copy()
    df['E_mean'] = E_mean
    df['E_std'] = E_std
    df['E_zscore'] = E_zscore
    df['f_direct'] = f_direct
    df['delta_f_used'] = delta_f
    
    # Redistribute remaining probability (Equations 8-9)
    remaining = 1.0 - f_direct
    total_indirect = p2_base + p3_base
    
    if total_indirect > 1e-10:
        p2_local = remaining * (p2_base / total_indirect)
        p3_local = remaining * (p3_base / total_indirect)
    else:
        # Carbon case: p2 ≈ 0, all indirect goes to p3
        p2_local = np.zeros_like(f_direct)
        p3_local = remaining
    
    df['p2_local'] = p2_local
    df['p3_local'] = p3_local
    
    # Reference p_ind at normoxia
    p_ind_ref = (21.0 + K_fix) / (21.0 + K_fix + K_repair)
    
    # Compute P_DSB for each O2 level (Equation 10)
    for O2 in O2_levels:
        p_ind = (O2 + K_fix) / (O2 + K_fix + K_repair)
        
        # Raw retention probability
        P_raw = f_direct + p2_local * p_ind + p3_local * p_ind**2
        P_raw_ref = f_direct + p2_local * p_ind_ref + p3_local * p_ind_ref**2
        
        # Normalized P_DSB (ensures P_DSB = 1 at normoxia)
        P_DSB = np.where(P_raw_ref > 0, P_raw / P_raw_ref, 1.0)
        P_DSB = np.clip(P_DSB, 0, 1)
        
        # Column naming: replace decimal point with 'p'
        col_name = f'P_DSB_{O2:.3f}'.replace('.', 'p')
        df[col_name] = P_DSB
    
    return df


# STATISTICS AND OUTPUT

def compute_summary_statistics(df: pd.DataFrame, particle_type: str, dataset_type: str = "calibration") -> Dict:
    """Compute comprehensive summary statistics."""
    stats_dict = {
        'particle': particle_type,
        'dataset_type': dataset_type,
        'n_dsbs': len(df),
        'extraction_timestamp': datetime.now().isoformat(),
        'E_local': {
            'mean': float(df['E_local'].mean()),
            'std': float(df['E_local'].std()),
            'min': float(df['E_local'].min()),
            'max': float(df['E_local'].max()),
            'median': float(df['E_local'].median()),
            'q25': float(df['E_local'].quantile(0.25)),
            'q75': float(df['E_local'].quantile(0.75)),
            'cv': float(df['E_local'].std() / df['E_local'].mean() * 100) if df['E_local'].mean() > 0 else 0,
            'n_zero': int((df['E_local'] == 0).sum()),
            'n_nonzero': int((df['E_local'] > 0).sum()),
        },
        'E_total_3x3x3': {
            'mean': float(df['E_total_3x3x3'].mean()),
            'std': float(df['E_total_3x3x3'].std()),
            'min': float(df['E_total_3x3x3'].min()),
            'max': float(df['E_total_3x3x3'].max()),
        },
        'grad_mag': {
            'mean': float(df['grad_mag'].mean()),
            'std': float(df['grad_mag'].std()),
            'max': float(df['grad_mag'].max()),
        },
        'core_frac': {
            'mean': float(df['core_frac'].mean()),
            'std': float(df['core_frac'].std()),
        },
        'n_active_3x3x3': {
            'mean': float(df['n_active_3x3x3'].mean()),
        },
        'is_local_max': {
            'count': int(df['is_local_max'].sum()),
            'fraction': float(df['is_local_max'].mean()),
        },
    }

    if 'complexity' in df.columns:
        stats_dict['complexity_distribution'] = {
            str(k): int(v) for k, v in df['complexity'].value_counts().items()
        }
    
    # Add P_DSB statistics if computed
    pdsb_cols = [c for c in df.columns if c.startswith('P_DSB_')]
    if pdsb_cols:
        stats_dict['P_DSB'] = {}
        for col in pdsb_cols:
            O2_str = col.replace('P_DSB_', '').replace('p', '.')
            stats_dict['P_DSB'][O2_str] = {
                'mean': float(df[col].mean()),
                'std': float(df[col].std()),
                'cv': float(df[col].std() / df[col].mean() * 100) if df[col].mean() > 0 else 0,
                'min': float(df[col].min()),
                'max': float(df[col].max()),
            }

    return stats_dict


def print_summary_table(all_stats: Dict[str, Dict]) -> None:
    """Print formatted summary table."""
    print("\n" + "=" * 95)
    print("ENERGY FEATURE EXTRACTION SUMMARY (VOxA Model)")
    print("=" * 95)

    print(f"\n{'Particle':<12} {'Dataset':<12} {'N DSBs':<10} {'E_local mean':<14} {'E_local CV':<12} "
          f"{'N E>0':<10}")
    print("-" * 95)

    for particle, stats in all_stats.items():
        e = stats['E_local']
        print(f"{particle:<12} {stats.get('dataset_type', 'unknown'):<12} {stats['n_dsbs']:<10} "
              f"{e['mean']:<14.6e} {e['cv']:<12.1f}% {e['n_nonzero']:<10}")

    print("=" * 95)


# MAIN PROCESSING

def process_particle_data(
    dsb_file: Path,
    energy_file: Path,
    particle_type: str,
    output_dir: Path,
    coord_columns: Tuple[str, str, str] = ('x_um', 'y_um', 'z_um'),
    geom: SimulationGeometry = None,
    run_id: str = "",
    dataset_type: str = "calibration",
    compute_predictions: bool = True
) -> Tuple[pd.DataFrame, Dict]:
    """Process a single particle's data files."""
    if geom is None:
        geom = SimulationGeometry()

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    logger.info(f"\n{'='*70}")
    logger.info(f"Processing {particle_type.upper()} ({dataset_type})")
    logger.info(f"{'='*70}")

    # Parse files
    dsb_df = parse_dsb_file(dsb_file, coord_columns)
    energy_grid, energy_meta = parse_topas_energy_file(energy_file)

    # Validate geometry
    if energy_meta.get('bin_size_cm'):
        geom.validate_against_topas_header(energy_meta['bin_size_cm'], 'energy')

    # Extract features
    logger.info("Extracting energy features for each DSB...")
    features = extract_dsb_energy_features(
        dsb_df, energy_grid, coord_columns,
        complexity_column='complexity',
        geom=geom,
        particle=particle_type,
        run_id=run_id,
        dataset_type=dataset_type
    )
    features_df = features_to_dataframe(features)

    # Compute VA predictions
    if compute_predictions:
        logger.info("Computing VOxA VA predictions...")
        features_df = compute_va_predictions(features_df, particle_type)

    # Compute statistics
    stats = compute_summary_statistics(features_df, particle_type, dataset_type)
    stats['energy_file_metadata'] = energy_meta

    # Save outputs
    prefix = f"{particle_type}_{dataset_type}"
    features_csv = output_dir / f"{prefix}_dsb_energy_features.csv"
    features_df.to_csv(features_csv, index=False)
    logger.info(f"Features saved to: {features_csv}")

    stats_json = output_dir / f"{prefix}_energy_statistics.json"
    with open(stats_json, 'w') as f:
        json.dump(stats, f, indent=2, default=str)
    logger.info(f"Statistics saved to: {stats_json}")

    return features_df, stats


# MAIN

def main():
    parser = argparse.ArgumentParser(
        description="Energy Feature Extraction for VOxA Model",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Modes:
  calibration  Extract features from 2500-DSB calibration datasets (Step 9 prep)
  validation   Extract features from 400-DSB validation datasets (Step 11 prep)

Examples:
  # Extract calibration features
  python extract_energy_features_voxa.py --mode calibration \\
      --data-dir ./calibration_data \\
      --output-dir ./voxa_features_output_calibration

  # Extract validation features
  python extract_energy_features_voxa.py --mode validation \\
      --data-dir ./validation_data \\
      --output-dir ./voxa_features_output_validation
        """
    )
    
    parser.add_argument('--mode', '-m', type=str, default='calibration',
                       choices=['calibration', 'validation'],
                       help='Processing mode: calibration (2500-DSB) or validation (400-DSB)')
    parser.add_argument('--data-dir', '-d', type=str, default='./data',
                       help='Directory containing DSB and energy files')
    parser.add_argument('--output-dir', '-o', type=str, required=True,
                       help='Output directory for extracted features (e.g., voxa_features_output_calibration)')
    parser.add_argument('--run-id', '-r', type=str, default='',
                       help='Run identifier for tracking')
    parser.add_argument('--no-predictions', action='store_true',
                       help='Skip VOxA P_DSB predictions')
    parser.add_argument('--particles', '-p', type=str, nargs='+',
                       default=['electron', 'proton', 'carbon'],
                       help='Particles to process')

    args = parser.parse_args()

    print("\n" + "=" * 70)
    print("ENERGY FEATURE EXTRACTION FOR VOxA MODEL")
    print(f"Mode: {args.mode.upper()}")
    print("=" * 70)
    
    # Feature extraction mode
    data_dir = Path(args.data_dir)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    print(f"Data directory: {data_dir}")
    print(f"Output directory: {output_dir}")
    print(f"Dataset type: {args.mode}")
    print(f"Particles: {args.particles}")
    print(f"Compute predictions: {not args.no_predictions}")

    # File patterns to try
    file_patterns = {
        'dsb': [
            '{particle}_dsb_complexity.csv',
            '{particle}_DSB_complexity.csv',
            '{particle}_dsb.csv',
            '{particle}_DSBs.csv',
        ],
        'energy': [
            '{particle}_EnergyDeposit.csv',
            '{particle}_energy_deposit.csv',
            '{particle}_Energy.csv',
        ]
    }

    # Find files for each particle
    particle_files = {}
    for particle in args.particles:
        dsb_file = None
        energy_file = None
        
        for pattern in file_patterns['dsb']:
            path = data_dir / pattern.format(particle=particle)
            if path.exists():
                dsb_file = path
                break
        
        for pattern in file_patterns['energy']:
            path = data_dir / pattern.format(particle=particle)
            if path.exists():
                energy_file = path
                break
        
        if dsb_file and energy_file:
            particle_files[particle] = {'dsb': dsb_file, 'energy': energy_file}
        else:
            if not dsb_file:
                logger.warning(f"DSB file not found for {particle}")
            if not energy_file:
                logger.warning(f"Energy file not found for {particle}")

    if not particle_files:
        print("\nNo valid particle data found!")
        print(f"Please ensure files exist in: {data_dir}")
        print("\nExpected file patterns:")
        for ftype, patterns in file_patterns.items():
            print(f"  {ftype}: {patterns[0]}")
        return

    # Process each particle
    geom = SimulationGeometry()
    geom.print_summary()

    coord_columns = ('x_um', 'y_um', 'z_um')
    
    all_features = {}
    all_stats = {}

    for particle, files in particle_files.items():
        try:
            features_df, stats = process_particle_data(
                dsb_file=files['dsb'],
                energy_file=files['energy'],
                particle_type=particle,
                output_dir=output_dir,
                coord_columns=coord_columns,
                geom=geom,
                run_id=args.run_id,
                dataset_type=args.mode,
                compute_predictions=not args.no_predictions
            )
            all_features[particle] = features_df
            all_stats[particle] = stats
            
        except Exception as e:
            logger.error(f"Error processing {particle}: {e}")
            import traceback
            traceback.print_exc()

    # Generate combined output
    if all_features:
        # Combined features CSV
        dfs = []
        for particle, df in all_features.items():
            dfs.append(df)
        
        combined = pd.concat(dfs, ignore_index=True)
        combined_path = output_dir / f"all_particles_{args.mode}_energy_features.csv"
        combined.to_csv(combined_path, index=False)
        logger.info(f"Combined features saved to: {combined_path}")
        logger.info(f"  Total DSBs: {len(combined)}")

        # Combined statistics JSON
        combined_stats = {
            'extraction_timestamp': datetime.now().isoformat(),
            'voxa_version': VOXA_REFERENCE['model_version'],
            'dataset_type': args.mode,
            'geometry': geom.to_dict(),
            'voxa_reference': VOXA_REFERENCE,
            'particles': all_stats,
            'total_dsbs': len(combined),
        }

        stats_path = output_dir / f"combined_{args.mode}_statistics.json"
        with open(stats_path, 'w') as f:
            json.dump(combined_stats, f, indent=2, default=str)
        logger.info(f"Combined statistics saved to: {stats_path}")

        print_summary_table(all_stats)

    # Final summary
    print("\n" + "=" * 70)
    print("EXTRACTION COMPLETE")
    print("=" * 70)
    print(f"\nOutput directory: {output_dir.absolute()}")
    print(f"Dataset type: {args.mode}")
    print(f"Total particles processed: {len(all_features)}")
    print(f"Total DSBs: {sum(len(df) for df in all_features.values())}")
    
    print("\nOutput files generated:")
    print(f"  • all_particles_{args.mode}_energy_features.csv")
    print(f"  • combined_{args.mode}_statistics.json")
    for particle in all_features.keys():
        print(f"  • {particle}_{args.mode}_dsb_energy_features.csv")
        print(f"  • {particle}_{args.mode}_energy_statistics.json")
    
    print("\nNext steps:")
    if args.mode == 'calibration':
        print("  1. Review combined_calibration_statistics.json")
        print("  2. Run Step 9: VA calibration")
        print("     Rscript step9_voxa_voxel_aware_calibration.R")
        print("  3. Run Step 10: VA diagnostics")
        print("     Rscript step10_voxa_va_diagnostics.R")
    elif args.mode == 'validation':
        print("  1. Review combined_validation_statistics.json")
        print("  2. Run Step 11: Scaling validation (requires both calibration and validation)")
        print("     Rscript step11_voxa_scaling_validation.R")
    print()


if __name__ == "__main__":
    main()