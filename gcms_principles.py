#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Test suite: Blind Deconvolution of 2D Signals using a Symbolic DNA Representation

Case Study: Identification of Non-Intentionally Added Substances (NIAS) in Recycled Polymers

This demonstration is part of the Generative Simulation Initiative
Contact: olivier.vitrac@gmail.com
"""

import os
import numpy as np
from sig2dna_core.signomics import signal_collection

# Setup for saving figures
outputfolder = "./images" if os.path.isdir("./images") else ("../images" if os.path.isdir("../images") else None)

# %% Step 1: Generate synthetic 2D GC-MS-like signals with controlled overlaps
#
# Each 2D signal is generated as a collection of m time-series (t,) that represent ion channels.
# The simulation mimics non-overlapping compounds whose retention (time) and fragmentation (m/z)
# properties are modeled using a sparse number of randomly placed peaks.
#
# Overlaps between substances are introduced through linear combinations.

t = 1024  # Number of time samples (dim1)
m = 32    # Number of ion channels or m/z values (dim2)
n_peaks = (6, 10)
n2Dsignals = 5  # Number of source signals (e.g., hypothetical pure substances)

for i in range(n2Dsignals):
    sig = signal_collection.generate_synthetic(
        n_signals=m,           # Each signal is a row of the (t, m) 2D matrix
        n_peaks=n_peaks,       # Random number of peaks
        kinds=("gauss",),      # Shape of the peaks
        width_range=(0.5, 3),  # Peak widths
        height_range=(1.0, 5.0),  # Peak heights
        x_range=(0, t-1),
        n_points=t,
        normalize=False,
        seed=40 + i * 10,
        name_prefix=f"G{i}"
    )[0]
    gcms = sig if i == 0 else gcms + sig  # Create overlapping signals

# Visualize each 1D signal in the 2D matrix
fig1 = gcms.plot()

# %% Step 2: Convert each time signal to a symbolic DNA representation
#
# Using second-derivative logic at scale=4, each (t,) signal is encoded into a sequence of symbols.
# Each character encodes a local pattern or curvature change: A, B, C, X, Y, Z, or blank (_).
# This dense encoding preserves signal length and local morphology.

scale = 4
dna_gcms = gcms._toDNA(scales=scale)  # Wrap into a DNAsignal_collection

# Display a heatmap of the first 100 symbols in all m signals
print("Alphabet used in DNA encoding:", dna_gcms.letters)
fig2 = dna_gcms.plot_letters()

# %% Step 3: Sinusoidal encoding of symbolic DNA positions
#
# We encode each character position into a sinusoidal latent space using the "transformer" trick:
# $$ PE_t(u,d) = \sin(u \cdot f_d) $$
# This preserves order and allows geometrical projection of symbolic dynamics.

d = 128  # Latent dimensionality (should exceed number of latent sources)
dna_gcms.sinencode_dna_full(d_model=d, operation="sum")

# Visualize letter-wise latent space embeddings (1D projections)
fig3 = dna_gcms.plot()

# Latent space projection of letter embeddings (dimension reduction view)
fig4 = dna_gcms.plot_embedding_projection()

# %% Step 4: Full 2D symbolic encoding with spatial semantics
#
# In the general model:
#     $$ v_{t,m,d} = E_{t,m,d} + PE_t(t,d) + PE_m(m,d) $$
# where PE_t and PE_m encode positions in time and m/z space, respectively.
# This version keeps 2D semantic structure for wrapped encoding (not flattened).

fig5 = dna_gcms.plot_v_symbol_components()

# %% Step 5: Raster scan simplification for GC-MS signals
#
# GC-MS signals are multiplexed in time, and ion channels are not acquired in parallel.
# The simplified raster-based encoding model removes PE_m:
#     $$ v_{u,d} = E_{u,d} + PE_t(u,d) $$
# where u is the linear index over the 2D matrix (raster scan).
# This variant reduces memory use and simplifies decoding.

fig6 = dna_gcms.plot_vtm_full()

# %% Step 6: Blind Deconvolution via PCA
#
# The goal is to decompose the encoded 2D signal into orthogonal components:
#     $$ X = \sum_k \sigma_k C_k $$
# where $C_k$ are PCA eigenvectors and $\sigma_k$ are projections onto latent dimensions.
# This identifies sources without prior knowledge of peak positions.

figsources = dna_gcms.deconvolve_latent_sources()[3]

# %% Step 7: Save figures
fig1.print("gcms1_Flattened_2Dsignals", outputfolder)
fig2.print("gcms2_Letters_heatmap", outputfolder)
fig3.print("gcms3_Letters_sinencoded", outputfolder)
fig4.print("gcms4_Letters_projections", outputfolder)
fig5.print("gcms5_Letters_full2D", outputfolder)
fig6.print("gcms6_Letters_raster2D", outputfolder)

for i, f in enumerate(figsources.keys()):
    figsources[f].print(f"gcms7{chr(i + 97)}_deconvolution_{f}", outputfolder)
