#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Test suite for the blind deconvolution of 2D signals using a DNA representation of signals.'d'

Case study: Indentication of non-intentionally added substances in recycled polymers

Generative Simulation Initiative\olivier.vitrac@gmail.com
"""

import os
import numpy as np
from sig2dna_core.signomics import signal_collection

outputfolder = "./images" if os.path.isdir("./images") else ("../images" if os.path.isdir("../images") else None)
# %% generate a linear combination of overlaping 2D random signals
# These 2D signals mimic GC-MS spectra but without chemical meaning.
# Indeed, the peaks are random in both directions (dim1=time, dim2=m or m/z).
#
# A 2D signal is a collection (m-list) of (t,) signals representing a (t,m) 2D signal.
# Linear operations + and * can be applied to compatible (t,m) collections
# To use the analogy with GC-MS:
#   gcms[i] represents the (i+1)th IC
#   m=len(gcms) is the mass spectrum resolution (number of m or m/z)
#
# The number of inserted peaks is controlled by n_peaks

t = 1024 # number of sampling times (dim1)
m = 32  # number of m/z or number of IC (dim2)
n_peaks = (6,10)
n2Dsignals = 5 # number of dependent 2D signals

for i in range(n2Dsignals):
    # non-overlapping 2D signals as a collection
    tmp = signal_collection.generate_synthetic(
            n_signals=m,             # number of (t,) signals (dim2)
            n_peaks=n_peaks,         # 8 will force eight peaks, (4,12) codes from 4 to 12 peaks
            kinds=("gauss",),
            width_range=(0.5, 3),    # range of peak withs
            height_range=(1.0, 5.0), # range of peak heights
            x_range=(0, t-1),        # x is the index (dim1)
            n_points=t,              # x is the index (dim1)
            normalize=False,
            seed=40+i*10,
            name_prefix=f"G{i}")[0]
    if i==0:
        gcms = tmp # initialization
    else:
        gcms += tmp # create overlaps by adding dependent signals (mathematical addition)

# control plots of all signals (flattened)
fig1 = gcms.plot()

# %% Conversion of analogic 2D signals into 2D a DNAsignal-collection (2D-wrapped DNAstr signals)
#
# This conversion is essentially the sign variation of the second derivative of the temporal (1D)
# original signal at scale=4. The numeric signal is converted into strings (DNAstr) consisting of
# an alphabet of 6 letters ('A', 'B', 'C', 'X', 'Y', 'Z') + blank ("_").
# The dense (full) representation of the code guarantees that all transformed (t,) signals keep the
# same length t. The position of each letter will be subsequently used to characterize each signal.

scale = 4 # scale for the continuous-walet-transform analysis
dna_gcms = gcms._toDNA(scales=scale) # convert GC-MS signals to DNAsignals

# Plot letters as a 2D heatmap
print("letters=",dna_gcms.letters) # letters
fig2 = dna_gcms.plot_letters() # their positions as a heatmap



# %% Sinusoidal encoding of the positions of letters in 2D-wrapped signals
#
# The transformation to DNAsignal at scale and its sinusoidal encoding assumes that
# signals are sequentually acquired (multiplexed) so that no specific sinusoidal encoding
# is required along dim2.
#
# dna_gcms.rasterscan = True is implicitely assumed, reducing the memory footprint.
#
# The positions of each letter ['A', 'B', 'C', 'X', 'Y', 'Z', '_'] are encoded independently in
# a sinusoidal space with d dimensions, where d should be chosen significantly greater than the
# number of latent sources to be present in the 2D mixture.

d = 128   # number of latent dimensions (sinusoidal encoding)
# sum is the reduction method to keep only transformed signals and their subspace
dna_gcms.sinencode_dna_full(d_model=d,operation="sum")


# Plot the encoded positions of each letter (reduced signals)
fig3 = dna_gcms.plot()         # shows 1D-letter signal positions encoded in sinusoidal space (latent space)

# Separate the m-DNAsignals based on the sole sinusoidal ecoding of each letter (i.e. reduced signals)
# Each (m,) signal is compressed to (d,) signal
fig4 = dna_gcms.plot_embedding_projection() # shows each letters in a 2D projection of the latent space

# Full 2D coding of positions using the general 2D model (no flattening assumption)
#$$
#    v_{u,m,d} = E_{u,m,d} + PE_t(u,d) + PE_m(m,d)
#$$
# where
#    - $u$ is the unique index of the 2D signal (u=0..t*m-1)
#    - $E_{u,m,d}$: symbolic encoding across segments.letters.
#    - $PE_t$: positional encoding along the time/segment axis (t).
#    - $PE_m$: positional encoding along the mass channel/ion axis (m).
fig5 = dna_gcms.plot_v_symbol_components()

# Note:
# The previous approach does not preserves the interpretability of transformed signals and the original
# meaning (e.g. GC-MS) cannot be restored. Additional assumptions and simplications are required to keep
# the meaningfull order of letters.

# Simplication for multiplexed GC-MS signals (all measurements are carried out sequentially).
# PE_m(m,d) has no physical meaning in this case (we do not need to encode the position in dim2), it
# can be dropped.
# The central idea is :
#   1) to drift/translate the letter position in the sinusoidal latent space by a (1,d) vector that
#      represents the letter/segment in the DNAstr subsignal.
#   2) the (1,d) vector representing the letter is chosen randomly not to emphasize any dimension in
#      the latent space.
#   3) As before, all positions refers to the origin of the (t,) subsignal (dim1)
#
# The revised model reads:
#$$
#    v_{u,d} = E_{u,d} + PE_t(u,d)
#$$
# where
#   - $u$ is as before the unique index of the 2D signal (u=0..t*m-1)
#   - $E_{u,d}$ represents the coding of all letters ("A","B","Y","Z"...) in their original order
#     as random(1,d) vectors.
#   - $PE_t(u,:)$: represents the sinusoidal encoding
fig6 = dna_gcms.plot_vtm_full()

# %% Blind deconvolution of the 2D transformed signal (e.g. separation of 2D chromatograms)
# Each v_{u,k} represents the chromatrogram at the given sinusoidal scale k=0..d-1 regardless of the
# intensity of the peaks. This redundant information can be used to recover chromatograms with
# non-overlapping peaks (along dim 1). This operation can be carried out with a PCA or a UMAP.
# The best analogy is to imagine that two texts with the same alphabet have been printed on the same
# page and we try to recover the idenpendently the two overlapping texts without understanding the
# language used to write the text.

figsources = dna_gcms.deconvolve_latent_sources()[3]


# %% printouts
fig1.print("gcms1_Flatened_2Dsignals", outputfolder)
fig2.print("gcms2_Letters_heatmap", outputfolder)
fig3.print("gcms3_Letters_sinencoded", outputfolder)
fig4.print("gcms4_Letters_projections", outputfolder)
fig5.print("gcms5_Letters_full2D", outputfolder)
fig6.print("gcms6_Letters_raster2D", outputfolder)
for i,f in enumerate(figsources.keys()):
    figsources[f].print(f"gcms7{chr(i+96)}_{f}_variance", outputfolder)
