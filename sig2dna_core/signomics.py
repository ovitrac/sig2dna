#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Module: sig2dna_core.signomics.py (from the Generative Simulation Initiative)
==================================

This is the core module of the sig2dna framework, dedicated to transforming numerical chemical signals into DNA-like symbolic representations. The module enables symbolic analysis, fingerprinting, alignment, and classification of complex analytical signals, such as:

- GC-MS / GC-FID
- HPLC-MS
- NMR / FTIR / Raman
- RX and other spectroscopy data

It is designed to facilitate high-throughput pattern recognition, compression, clustering, and AI/ML-based classification. Symbolic transformation is based on **wavelet decomposition (Mexican Hat/Ricker)** and segment encoding into letters (e.g., A, B, C, X, Y, Z). The representation preserves key structural patterns (e.g., peak transitions) across multiple scales and supports entropy-based distances.

Main Components
---------------
- DNAsignal — core class to transform a signal into DNA-like symbolic representation
- DNAstr — string subclass enabling alignment, entropy analysis, visualization, and reconstruction
- DNApairwiseAnalysis — distance and clustering toolbox for aligned DNA codes (PCoA, dendrogram, 2D/3D plots)

Key Features
------------
- Multi-scale wavelet transform with symbolic encoding
- Symbolic entropy and mutual information measures
- Fast symbolic alignment using `difflib` or `biopython`
- Pairwise symbolic distances (Shannon, excess entropy, Jaccard, Jensen-Shannon)
- Interactive plotting: segments, alignment masks, triangle patches
- Motif search, alignment visualization, HTML and terminal rendering
- Dimensionality reduction (MDS), clustering, dendrogram and heatmaps

Core Concept
------------

Input Signal
- 1D NumPy array S of shape (m,)
- Data type: np.float64 (default) or np.float32
- Typically sparse and non-negative, such as GC-MS total ion chromatograms

Wavelet Transform
- CWT with Ricker (Mexican hat) wavelet
- Scales: $s = 2^0, 2^1,..., 2^n
- Downsampling by scale (to reduce data volume and capture features at relevant resolutions)

Symbolic Encoding and compressed representation

The transformed signal Ts at scale s is converted to a sequence of symbolic letters using the rules:

+--------+---------------------------------------------------------------+
| Symbol | Description                                                   |
+========+===============================================================+
| A      | Monotonic increase crossing from − to +                       |
+--------+---------------------------------------------------------------+
| B      | Monotonic increase from − to − (no zero crossing)             |
+--------+---------------------------------------------------------------+
| C      | Monotonic increase from + to +                                |
+--------+---------------------------------------------------------------+
| X      | Monotonic decrease from + to + (no zero crossing)             |
+--------+---------------------------------------------------------------+
| Y      | Monotonic decrease from − to −                                |
+--------+---------------------------------------------------------------+
| Z      | Monotonic decrease crossing from + to −                       |
+--------+---------------------------------------------------------------+
| _      | Zero or noise (after filtering)                               |
+--------+---------------------------------------------------------------+


Each encoded segment is associated with:
- width: number of points
- height: amplitude difference

These form the compressed representation

Installation
------------
Install all dependencies with:

conda install pywavelets seaborn scikit-learn
conda install -c conda-forge python-Levenshtein biopython

Examples
--------
    >>> from signomics import DNAsignal
    >>> from signal import signal

    >>> # Load a sampled signal (e.g., from GC-MS, Raman)
    >>> S = signal.from_peaks(...)  # or any constructor for sampled signals

    >>> # Encode into DNA-like format
    >>> D = DNAsignal(S, encode=True)
    >>> D.encode_dna()
    >>> D.plot_codes(scale=4)

    >>> # Compare samples and cluster
    >>> Dlist = [DNAsignal(S1, encode=True), DNAsignal(S2, encode=True), ...]
    >>> analysis = DNAsignal._pairwiseEntropyDistance(Dlist, scale=4)
    >>> analysis.plot_dendrogram()
    >>> analysis.scatter(n_clusters=3)

Notes
-----
The methodology implemented in this module covers and extends the approaches initially tested during the PhD of Julien Kermorvant.
"Concept of chemical fingerprints applied to the management of chemical risk of materials, recycled deposits and food packaging". PhD thesis AgroParisTech. December 2023. https://theses.hal.science/tel-04194172


Maintenance & forking
---------------------

        $ git init -b main
        $ gh repo create sig2dna --public --source=. --remote=origin --push
        $ # alternatively
        $ # git remote add origin git@github.com:ovitrac/sig2dna.git
        $ # git branch -M main        # Ensure current branch is named 'main'
        $ # git push -u origin main   # Push and set upstream tracking

        $ tree -P '*.py' -P '*.md' -P 'LICENSE' -I '__pycache__|.*' --prune
        $ conda activate base
        $ pdoc ./sig2dna_core/signomics.py -f --html -o ./docs
        $ doctoc --github  --maxlevel 2 README.md

        $ conda activate sphinxdoc
        $ cd docs_sphinx/
        $ make clean
        $ make html
        $ cp -rp build/html/. ../docs


Author: Olivier Vitrac — olivier.vitrac@gmail.com
Revision: 2025-06-13
"""

# %% Indentication
__project__ = "Signomics"
__author__ = "Olivier Vitrac"
__copyright__ = "Copyright 2025"
__credits__ = ["Olivier Vitrac"]
__license__ = "MIT"
__maintainer__ = "Olivier Vitrac"
__email__ = "olivier.vitrac@gmail.com"
__version__ = "0.51"

# %% Dependencies
# note:
# cwt and ricker are depreciated since scipy v1.12, pwywalets is used instead. See: https://docs.scipy.org/doc/scipy-1.12.0/reference/generated/scipy.signal.cwt.html
#
# Conda users (from a fresn environment)
# conda install -c conda-forge tqdm numpy pandas scipy matplotlib ipython ipykernel pywavelets python-Levenshtein biopython scikit-learn, seaborn

# Generic libs
import os, sys, socket, getpass, datetime, uuid, operator, json, gzip, hashlib, re, math
import inspect, importlib.util, warnings
from pathlib import Path
from collections import defaultdict, UserDict, OrderedDict, Counter
from collections.abc import Sequence
from types import SimpleNamespace
#from dataclasses import dataclass
#from itertools import islice
from copy import deepcopy
from time import time
from tqdm import tqdm

# Math, machine-learning and visualization libs
import random, pickle
import numpy as np
import pandas as pd
from numpy.lib.stride_tricks import sliding_window_view
from scipy.interpolate import interp1d
from scipy.signal import  medfilt
from scipy.ndimage import uniform_filter1d
from scipy.stats import entropy
from scipy.spatial.distance import jensenshannon
from scipy.optimize import minimize_scalar
from scipy.linalg import svd, pinv
from difflib import SequenceMatcher
import matplotlib.pyplot as plt
import matplotlib.cm as cm # from matplotlib.cm import get_cmap
import matplotlib.colors as mcolors
from matplotlib.patches import Polygon, Rectangle
from matplotlib.colors import Normalize, ListedColormap
import matplotlib.patches as mpatches
from mpl_toolkits.mplot3d import Axes3D
from IPython.display import HTML, display

# Non-default libs
# Check availability of optional packages
optional_dependencies = {
    "pywt": {"package": "pywavelets", "optional": False, "source":""},
    "Levenshtein": {"package": "python-Levenshtein", "optional": False, "source":"-c conda-forge"},
    "Bio.Align": {"package": "biopython", "optional": True, "source": "-c conda-forge"},
    "seaborn": {"package": "seaborn", "optional": False, "source":""},
    "sklearn": {"package": "scikit-learn", "optional":False, "source":""},
    "umap":{"package": "umap-learn", "optional":True, "source":"-c conda-forge"}
}
dependency_status = []
for module, meta in optional_dependencies.items():
    module_name = module.split('.')[0]
    available = importlib.util.find_spec(module_name) is not None
    status = {
        "Module": module,
        "Available": available,
        "Install Instruction": f"conda install {meta['source']} {meta['package']}"
    }
    dependency_status.append(status)
    if not available:
        msg = (
            f"Required module '{module}' not found.\n"
            f"To install: {status['Install Instruction']}"
        )
        if meta.get("optional", False):
            warnings.warn(msg)
        else:
            raise ImportError(msg)

import pywt
import Levenshtein
from Bio.Align import PairwiseAligner
from sklearn.decomposition import PCA
from sklearn.manifold import MDS
from sklearn.metrics import pairwise_distances, silhouette_score
from sklearn.preprocessing import StandardScaler
from scipy.cluster.hierarchy import linkage, dendrogram, fcluster
from scipy.spatial.distance import squareform
import seaborn as sns
#from umap import UMAP

__all__ = ['DNACodes', 'DNAFullCodes', 'DNApairwiseAnalysis', 'DNAsignal', 'DNAsignal_collection', 'DNAstr', 'SinusoidalEncoder', 'generator', 'import_local_module', 'peaks', 'signal', 'signal_collection']

# %% load dynamically figprint without pythonpath modification or installation
def import_local_module(name: str, relative_path: str):
    """
    Import a module by name from a file path relative to the calling module.

    Usage here:
    -----------
    import_local_module("figprint", "figprint.py") replaces import figprint (zero-installation)

    Parameters:
    ----------
    name : str
        Name to assign to the module (used internally).
    relative_path : str
        Path relative to the calling module's location.

    Returns:
    -------
    module
        Imported module object.

    Example:
    --------
    >>> figprint = import_local_module("figprint", "figprint.py")
    >>> figprint.print_pdf(...)
    """
    # Get the caller's file path
    caller_frame = inspect.stack()[1]
    caller_path = os.path.abspath(os.path.dirname(caller_frame.filename))
    # Resolve full path to the module
    full_path = os.path.abspath(os.path.join(caller_path, relative_path))
    # Import as custom name
    spec = importlib.util.spec_from_file_location(name, full_path)
    if spec is None or spec.loader is None:
        raise ImportError(f"Cannot load module {name} from {full_path}")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod

# use figprint = import_local_module if you want to access directly print_pdf, print_png...
import_local_module("figprint", "figprint.py")


# %% Specialized Encoder
class SinusoidalEncoder:
    """
    🌀 Generic sinusoidal encoder/decoder supporting symbolic and numeric sequences.

    Each scalar value is transformed into a vector of dimension `d_model`, where alternating components
    contain sinusoidal features of increasing frequency. The mapping is based on:

        For k = 0 to d_model/2 - 1:
            f_{2k}(x)   = sin(x / r_k)
            f_{2k+1}(x) = cos(x / r_k)

        where r_k = N^(2k / d_model)

    This representation preserves relative positions and scaling in a smooth, topologically faithful
    embedding space. The class supports multiple decoding strategies, scaling logic, residual control,
    and round-trip verification.

    Parameters
    ----------
    d_model : int
        Dimensionality of each sinusoidal encoding (must be even).
    N : int
        Frequency base for the positional encoding.
    dtype : np.dtype
        Output data type (default: np.float32).

    Attributes
    ----------
    d_model : int
        Embedding dimensionality.
    N : int
        Frequency base.
    dtype : np.dtype
        Data type for encoded output.
    _last_input_type : type
        Last input type passed to `encode`.
    _last_input_length : int
        Last input length passed to `encode`.
    _scale : float or None
        Scaling factor applied to normalize input values.
    _auto_scale_enabled : bool
        Whether autoscaling is enabled.
    _decode_residual_tolerance : float
        Tolerance for residual error checking in decode verification.

    Methods
    -------
    encode(values, scale=None)
        Encodes a sequence of values (scalar or symbolic) into sinusoidal embeddings.
    decode(embedding, method='least_squares', return_error=False)
        Decodes the embedding to the original values using the selected inverse method.
    fit_encoder(values, target_range=10.0)
        Automatically estimates and stores a scaling factor to normalize input values.
    set_decode_tolerance(tol)
        Sets the residual error threshold above which decoding results will raise a warning.
    verify_roundtrip(values, method='least_squares', scale='auto', verbose=True, return_details=False)
        Checks round-trip accuracy of encoding and decoding. Warns if residuals exceed tolerance.

    sinencode_dna_grouped(code, d_part, N)
        Encodes a `codes` entry (triplet-based segments grouped by letter).
    sindecode_dna_grouped(grouped, reference_code, d_part, N)
        Reconstructs code dictionary from grouped sinusoidal embeddings.
    sinencode_dnafull_grouped(dnafull, d_model, N)
        Encodes `codesfull` entries (strings) into grouped embeddings.
    sindecode_dnafull_grouped(grouped)
        Decodes grouped full embeddings back to `DNAstr`.

    Static Methods
    --------------
    to_complex(emb)
        Convert a sinusoidal embedding (sin, cos) into complex numbers using Euler's identity.
    complex_distance(emb1, emb2, norm='L2')
        Compute pointwise distances between two embeddings in the complex sinusoidal space.
    angle_difference(emb)
        Compute angular differences Δθ between consecutive elements of a sinusoidal embedding.
    phase_alignment(emb, ref)
        Align the phase of an embedding `emb` to a reference embedding `ref` using complex phase factors.
    pairwise_similarity(emb, metric='cosine')
        Compute a pairwise similarity or distance matrix ('cosine' or 'L2') between all elements.
    group_centroid(emb, labels=None, return_std=False)
        Compute the centroid (and optionally standard deviation) of groups in complex embedding space.
    phase_unwrap(emb, normalize=False)
        Perform phase unwrapping (à la Fourier) on the sinusoidal embedding, optionally normalized to [0, 1].

    Example (without scaling)
    --------------------------
    >>> s = SinusoidalEncoder(8, 100)  # poor encoder (8 dimensions, high N)
    >>> a = s.encode([0, 1, 1, 2, 2, 3, 4, 5, 6, 6])
    >>> s.decode(a)

    Output:
        [0.0,
         1.0000000072927564,
         1.0000000072927564,
         1.9999999989612762,
         1.9999999989612762,
         2.9999999881593866,
         3.9999997833243714,
         5.000000005427378,
         5.999999552487366,
         5.999999552487366]

    Notes:
        - Use lower N (e.g., N = 1000) to compress phase variation and allow larger input range.
        - Use scaling (via `fit_encoder()` or `scale=`) for large or high-resolution inputs.

    Example (with scaling)
    ----------------------
    >>> s2 = SinusoidalEncoder(128, 10000)
    >>> a2 = s2.encode([0, 1, 1, 2, 2, 3, 4, 5, 6, 6, 7, 7, 7, 8, 8, 8, 8, 16, 99, 130], scale=100)
    >>> s2.decode(a2)

    Output:
        [0.0,
         1.0000000127956146,
         1.0000000127956146,
         1.999999922004248,
         1.999999922004248,
         2.999999960015871,
         3.9999999840341935,
         5.000000076901334,
         5.999999935247362,
         5.999999935247362,
         6.999999909109446,
         6.999999909109446,
         6.999999909109446,
         8.000000084954152,
         8.000000084954152,
         8.000000084954152,
         8.000000084954152,
         16.000000431016925,
         99.00000148970207,
         129.99999823222058]

    Advanced Example
    -----------------
    import numpy as np
    import matplotlib.pyplot as plt

    >>> # 1. Construct a test input signal with smooth and jump segments
    >>> x_smooth = np.linspace(0, 20, 100)
    >>> x_jumps = np.array([25, 25, 26, 27, 100, 101, 130])
    >>> x = np.concatenate([x_smooth, x_jumps])

    >>> # 2. Initialize encoder with high d_model and N
    >>> s = SinusoidalEncoder(d_model=128, N=10000)

    >>> # 3. Fit auto-scaling to compress input into sinusoidal-friendly space
    >>> s.fit_encoder(x, target_range=10)

    >>> # 4. Encode and decode using all robust methods
    >>> a = s.encode(x)
    >>> decoded_lsq, err_lsq = s.decode(a, method='least_squares', return_error=True)
    >>> decoded_svd, err_svd = s.decode(a, method='svd', return_error=True)

    >>> # 5. Compare errors
    >>> true = x
    >>> lsq_error = np.abs(decoded_lsq - true)
    >>> svd_error = np.abs(decoded_svd - true)

    >>> # 6. Plot results
    >>> fig, axs = plt.subplots(2, 2, figsize=(12, 8))
    >>> axs[0, 0].plot(true, label="Original")
    >>> axs[0, 0].plot(decoded_lsq, '--', label="Decoded (LSQ)")
    >>> axs[0, 0].plot(decoded_svd, ':', label="Decoded (SVD)")
    >>> axs[0, 0].set_title("Decoded vs Original")
    >>> axs[0, 0].legend()

    >>> axs[0, 1].plot(lsq_error, label="Abs Error (LSQ)")
    >>> axs[0, 1].plot(svd_error, label="Abs Error (SVD)")
    >>> axs[0, 1].set_yscale('log')
    >>> axs[0, 1].set_title("Absolute Decoding Error (log scale)")
    >>> axs[0, 1].legend()

    >>> axs[1, 0].plot(err_lsq, label="Residual Norm (LSQ)")
    >>> axs[1, 0].plot(err_svd, label="Residual Norm (SVD)")
    >>> axs[1, 0].set_yscale('log')
    >>> axs[1, 0].set_title("Reconstruction Residuals")
    >>> axs[1, 0].legend()

    >>> axs[1, 1].hist(lsq_error, bins=50, alpha=0.7, label="LSQ")
    >>> axs[1, 1].hist(svd_error, bins=50, alpha=0.5, label="SVD")
    >>> axs[1, 1].set_title("Histogram of Absolute Errors")
    >>> axs[1, 1].legend()

    >>> plt.suptitle("🌀 SinusoidalEncoder: Accuracy Evaluation", fontsize=14)
    >>> plt.tight_layout()
    >>> plt.show()

    References (for the encoding)
    ------------------------------
    * Vaswani et al. (2017), "Attention is All You Need"
    * https://en.wikipedia.org/wiki/Transformer_(deep_learning_architecture)
    """

    def __init__(self, d_model=96, N=10000, dtype=np.float32):
        """
        SinuosidalEncoder Constructor

        Parameters
        ----------
        d_model : int
            Dimensionality of each sinusoidal encoding (must be even).
        N : int
            Frequency base for positional encoding.
        dtype : np.dtype
            Output data type (default: np.float32).
        """
        if d_model % 2 != 0:
            raise ValueError("d_model must be even")
        self.d_model = d_model
        self.N = N
        self.dtype = dtype
        self._last_input_type = None
        self._scale = None  # auto-scaling factor
        self._auto_scale_enabled = False
        self._decode_residual_tolerance = 1e-3

    def set_decode_tolerance(self, tol=1e-3):
        """
        Set maximum acceptable residual error for decoding.

        Parameters
        ----------
        tol : float
            Residual threshold above which a warning is triggered.
        """
        self._decode_residual_tolerance = tol

    def fit_encoder(self, values, target_range=10.0):
        """
        Fit a scaling factor to normalize values into a target sinusoidal-safe range.

        Parameters
        ----------
        values : array-like
            Original values to encode (will determine scale).
        target_range : float
            Maximum scaled range to span (e.g. [0, 10]).

        Returns
        -------
        float
            Recommended scale factor stored internally.
        """
        values = np.asarray(values, dtype=float).flatten()
        span = np.max(values) - np.min(values)
        if span == 0:
            self._scale = 1.0
        else:
            self._scale = span / target_range
        self._auto_scale_enabled = True
        return self._scale

    def encode(self, values, scale=None):
        """
        Encode input values into sinusoidal embeddings.

        Parameters
        ----------
        values : array-like
            Values to encode.
        scale : float or None
            Rescaling factor. If None, auto-scaling is applied if enabled.

        Returns
        -------
        np.ndarray
            Embedded values of shape (n, d_model)
        """
        arr = self._parse(values)
        self._last_input_type = type(values)
        self._last_input_length = len(values)

        if scale is not None:
            self._scale = scale
        elif self._auto_scale_enabled and self._scale is not None:
            pass  # use previously fitted scale
        else:
            self._scale = None

        if self._scale:
            arr = arr / self._scale

        return self._sin_embed(arr)


    def decode(self, embedding, method='least_squares', return_error=False):
        """
        Decode sinusoidal embeddings back to original values using selected method.

        Parameters
        ----------
        embedding : np.ndarray
            Encoded sinusoidal array of shape (n, d_model)
        method : str
            Decoding strategy: 'least_squares' (default), 'optimize', or 'naive'
        return_error : bool
            If True, returns (decoded_values, residual_error) as a tuple

        Returns
        -------
        decoded : list or array
            Reconstructed input values
        residual : np.ndarray, optional
            Residuals of decoding (per sample), returned only if return_error=True
        """
        if self._last_input_type is None:
            raise RuntimeError("Cannot decode: no prior encoding available.")
        if method == 'least_squares' or method == "LS":
            values, error = self._decode_lsq(embedding)
        elif method == 'svd':
            values, error = self._decode_svd(embedding)
        elif method == 'optimize':
            values, error = self._decode_optimize(embedding)
        elif method == 'naive':
            values = self._decode_naive(embedding)
            error = np.full_like(values, np.nan)
        else:
            raise ValueError("Unknown method. Use 'least_squares', 'optimize', or 'naive'.")
        if self._scale is not None:
            values *= self._scale
        if return_error and hasattr(error, "__len__"):
            bad = np.where(error > self._decode_residual_tolerance)[0]
            if len(bad):
                print(f"⚠️  Warning: {len(bad)} / {len(error)} decoded values exceed residual threshold "
                      f"({self._decode_residual_tolerance}).")
        reconstructed = self._reconstruct_type(values)
        return (reconstructed, error) if return_error else reconstructed

    def _parse(self, values):
        """Parse values"""
        if isinstance(values, str):
            return np.arange(len(values))
        elif isinstance(values, DNAstr):
            return np.arange(len(values))
        elif isinstance(values, (list, tuple, np.ndarray)):
            return np.array(values, dtype=float).flatten()
        else:
            raise TypeError(f"Unsupported input type: {type(values)}")

    def _reconstruct_type(self, arr):
        """Reconstruct types"""
        if self._last_input_type is str:
            return ''.join(chr(65 + int(round(x)) % 26) for x in arr)
        elif self._last_input_type is list:
            return arr.tolist()
        elif self._last_input_type is tuple:
            return tuple(arr.tolist())
        elif self._last_input_type is np.ndarray:
            return arr.astype(float)
        else:
            raise TypeError("Cannot decode unknown or unsupported type")

    def _sin_embed(self, values):
        """
        Project scalar values into sinusoidal embedding space.

        Parameters
        ----------
        values : np.ndarray
            1D array of scalar values to encode.

        Returns
        -------
        np.ndarray
            Array of shape (n, d_model) with sinusoidal features.
        """
        values = np.asarray(values, dtype=float).reshape(-1, 1)
        d_half = self.d_model // 2
        k = np.arange(d_half).reshape(1, -1)
        r = self.N ** (2 * k / self.d_model)
        angles = values / r
        return np.concatenate([np.sin(angles), np.cos(angles)], axis=1).astype(self.dtype)

    def _inverse_sin_embed(self, emb):
        """
        Inverse sinusoidal encoding (approximate).

        Parameters
        ----------
        emb : np.ndarray
            Encoded vectors (n, d_part*2)
        d_part : int
            Dimension of each embedding part
        N : int
            Frequency base

        Returns
        -------
        np.ndarray
            Approximate decoded scalar values
        """
        d_part = emb.shape[1] // 2
        angles = np.arctan2(emb[:, :d_part], emb[:, d_part:])
        k = np.arange(d_part).reshape(1, -1)
        r = self.N ** (2 * k / self.d_model)
        return np.mean(angles * r, axis=1)


    def _decode_svd(self, emb):
        """
        SVD-based phase unwrapping decoder with regularized least-squares.

        Parameters
        ----------
        emb : np.ndarray
            Encoded array of shape (n, d_model)

        Returns
        -------
        values : np.ndarray
            Decoded values
        residuals : np.ndarray
            Reconstruction errors (L2 norm)
        """
        d_half = self.d_model // 2
        k = np.arange(d_half)
        r = self.N ** (2 * k / self.d_model)
        freqs = 1 / r  # shape (d_half,)

        z = self.to_complex(emb)  # shape (n, d_half)
        theta = np.unwrap(np.angle(z), axis=0)  # shape (n, d_half)

        # Use pseudoinverse for robust fit
        F = freqs.reshape(1, -1)
        U, S, Vt = np.linalg.svd(F, full_matrices=False)
        pinv = (Vt.T @ np.diag(1 / S) @ U.T).reshape(-1)  # shape (d_half,)

        values = theta @ pinv
        fitted = np.outer(values, freqs)
        residual = np.linalg.norm(theta - fitted, axis=1)
        return values, residual


    def _decode_lsq(self, emb):
        """
        Least-squares phase unwrapping decoder (robust fast decoder).

        Parameters
        ----------
        emb : np.ndarray
            Encoded array of shape (n, d_model)

        Returns
        -------
        values : np.ndarray
            Decoded values
        residuals : np.ndarray
            Reconstruction errors (L2 norm)
        """
        d_half = self.d_model // 2
        k = np.arange(d_half)
        r = self.N ** (2 * k / self.d_model)
        freqs = 1 / r

        z = self.to_complex(emb)
        theta = np.unwrap(np.angle(z), axis=0)

        values = np.sum(theta * freqs, axis=1) / np.sum(freqs ** 2)
        # reconstruction error (projected back then subtract)
        fitted = np.outer(values, freqs)
        residual = np.linalg.norm(theta - fitted, axis=1)
        return values, residual

    def _decode_optimize(self, emb):
        """
        Decode each embedding vector via scalar minimization.

        Parameters
        ----------
        emb : np.ndarray
            Encoded array of shape (n, d_model)

        Returns
        -------
        values : np.ndarray
            Decoded values
        residuals : np.ndarray
            L2 residuals of the solution
        """

        d_half = self.d_model // 2
        k = np.arange(d_half)
        r = self.N ** (2 * k / self.d_model)

        def f_embed(x):
            x = np.asarray(x).reshape(-1, 1)
            angles = x / r.reshape(1, -1)
            return np.hstack([np.sin(angles), np.cos(angles)])

        def decode_single(e):
            def loss(x):
                return np.linalg.norm(f_embed(x) - e.reshape(1, -1))
            res = minimize_scalar(loss, bounds=(-100, 100), method='bounded')
            return res.x if res.success else np.nan, res.fun if res.success else np.inf

        results = [decode_single(e) for e in emb]
        values, errors = zip(*results)
        return np.array(values), np.array(errors)

    def _decode_naive(self, emb):
        """
        Simple decoder based on mean projection (coarse and periodic).

        Parameters
        ----------
        emb : np.ndarray
            Encoded sinusoidal array (n, d_model)

        Returns
        -------
        np.ndarray
            Coarse estimate of decoded values
        """
        d_half = self.d_model // 2
        angles = np.arctan2(emb[:, :d_half], emb[:, d_half:])
        k = np.arange(d_half).reshape(1, -1)
        r = self.N ** (2 * k / self.d_model)
        return np.mean(angles * r, axis=1)

    @staticmethod
    def sinencode_dna_grouped(code, d_part=32, N=10000):
        """
        Encode symbolic code segments grouped by letter into sinusoidal embeddings.

        Parameters
        ----------
        code : dict
            Must contain:
                - 'letters': str
                - 'xloc': list of (start, end) tuples
                - 'widths': list of float
                - 'heights': list of float
        d_part : int
            Number of dimensions per field (position, width, height)
        N : int
            Sinusoidal frequency base

        Returns
        -------
        dict
            Dictionary of embeddings by letter: {letter: np.ndarray(n, 3*d_part)}
        """
        enc = SinusoidalEncoder(d_model=d_part, N=N)
        grouped = {}
        dx = code.get("dx", 1.0)
        for letter in set(code["letters"]):
            indices = [i for i, l in enumerate(code["letters"]) if l == letter]
            starts = np.array([code["xloc"][i][0] for i in indices])
            widths = np.array([code["widths"][i] for i in indices])
            heights = np.array([code["heights"][i] for i in indices])
            emb_start = enc._sin_embed(starts)
            emb_width = enc._sin_embed(widths)
            emb_height = enc._sin_embed(heights)
            grouped[letter] = np.hstack([emb_start, emb_width, emb_height])
        return grouped

    @staticmethod
    def sindecode_dna_grouped(grouped, reference_code, d_part=32, N=10000):
        """
        Decode sinusoidal embeddings grouped by letter into symbolic segments.

        Parameters
        ----------
        grouped : dict
            Dictionary of letter: embeddings
        reference_code : dict
            Must include 'dx' (sampling resolution).
        d_part : int
            Number of dimensions per part (start, width, height)
        N : int
            Frequency base

        Returns
        -------
        dict
            Dictionary with keys: letters, widths, heights, xloc, iloc, dx
        """
        dx = reference_code["dx"]
        all_letters, all_xloc, all_widths, all_heights, all_iloc = [], [], [], [], []
        for letter, emb in grouped.items():
            n = emb.shape[0]
            x = SinusoidalEncoder._inverse_sin_embed(emb[:, :d_part], d_part, N)
            w = SinusoidalEncoder._inverse_sin_embed(emb[:, d_part:2*d_part], d_part, N)
            h = SinusoidalEncoder._inverse_sin_embed(emb[:, 2*d_part:], d_part, N)
            all_letters.extend([letter] * n)
            all_widths.extend(w)
            all_heights.extend(h)
            all_xloc.extend([(xi, xi + wi) for xi, wi in zip(x, w)])
            all_iloc.extend([(int(round(xi / dx)), int(round((xi + wi) / dx))) for xi, wi in zip(x, w)])
        return {
            "letters": ''.join(all_letters),
            "widths": all_widths,
            "heights": all_heights,
            "xloc": all_xloc,
            "iloc": all_iloc,
            "dx": dx
        }

    @staticmethod
    def to_complex(emb):
        """
        Convert sinusoidal embedding into a complex array using Euler's identity.

        Parameters
        ----------
        emb : np.ndarray
            Array of shape (..., 2d) where sin/cos pairs are stored.

        Returns
        -------
        np.ndarray
            Complex array of shape (..., d) with values exp(i * theta)
        """
        d = emb.shape[-1] // 2
        sin, cos = emb[..., :d], emb[..., d:]
        return cos + 1j * sin

    @staticmethod
    def complex_distance(emb1, emb2, norm='L2'):
        """
        Compute distance between two encoded arrays using complex projection.

        Parameters
        ----------
        emb1, emb2 : np.ndarray
            Encoded arrays of shape (..., 2d) to compare.
        norm : str
            'L2' for Euclidean norm, 'cos' for cosine angle distance.

        Returns
        -------
        np.ndarray
            Distance values per sample (1D array)
        """
        z1 = SinusoidalEncoder.to_complex(emb1)
        z2 = SinusoidalEncoder.to_complex(emb2)
        if norm == 'L2':
            return np.linalg.norm(z1 - z2, axis=-1)
        elif norm == 'cos':
            dot = np.real(np.sum(z1 * np.conj(z2), axis=-1))
            norm1 = np.linalg.norm(z1, axis=-1)
            norm2 = np.linalg.norm(z2, axis=-1)
            return 1 - (dot / (norm1 * norm2 + 1e-9))
        else:
            raise ValueError("Unknown norm type. Use 'L2' or 'cos'.")

    @staticmethod
    def angle_difference(emb):
        """
        Compute angular differences (∆θ) between consecutive embeddings.

        Parameters
        ----------
        emb : np.ndarray
            Encoded array of shape (n, 2d)

        Returns
        -------
        np.ndarray
            Array of shape (n-1, d) of angular differences in radians
        """
        z = SinusoidalEncoder.to_complex(emb)
        return np.angle(z[1:] * np.conj(z[:-1]))

    @staticmethod
    def phase_alignment(emb, ref):
        """
        Align embedding `emb` to reference `ref` using complex phase.

        Parameters
        ----------
        emb : np.ndarray
            Encoded array to align (n, 2d)
        ref : np.ndarray
            Reference encoded array (n, 2d)

        Returns
        -------
        np.ndarray
            Aligned encoding of `emb`, same shape as input.
        """
        z1 = SinusoidalEncoder.to_complex(emb)
        z2 = SinusoidalEncoder.to_complex(ref)
        phase_shift = np.sum(z2 * np.conj(z1), axis=-1, keepdims=True)
        align_factor = phase_shift / np.abs(phase_shift + 1e-9)
        z_aligned = z1 * align_factor
        return np.concatenate([np.imag(z_aligned), np.real(z_aligned)], axis=-1)

    @staticmethod
    def pairwise_similarity(emb, metric='cosine'):
        """
        Compute a pairwise similarity (or distance) matrix in sinusoidal embedding space.

        Parameters
        ----------
        emb : np.ndarray
            Encoded array of shape (n, 2d).
        metric : str
            Distance metric: 'cosine' for 1 - cosine similarity, 'L2' for Euclidean norm.

        Returns
        -------
        np.ndarray
            Pairwise similarity matrix of shape (n, n)
        """
        z = SinusoidalEncoder.to_complex(emb)
        if metric == 'cosine':
            norm = np.linalg.norm(z, axis=1, keepdims=True)
            Z = z / (norm + 1e-9)
            return 1 - np.real(Z @ Z.conj().T)
        elif metric == 'L2':
            diff = np.expand_dims(z, 1) - np.expand_dims(z, 0)
            return np.linalg.norm(diff, axis=-1)
        else:
            raise ValueError("Unknown metric: use 'cosine' or 'L2'")

    @staticmethod
    def group_centroid(emb, labels=None, return_std=False):
        """
        Compute the centroid (average embedding) of each group in complex sinusoidal space.

        Parameters
        ----------
        emb : np.ndarray
            Encoded array of shape (n, 2d), with sin and cos interleaved.
        labels : list or np.ndarray, optional
            Group labels (n,). If None, the entire set is treated as one group.
        return_std : bool
            Whether to also return the standard deviation per group.

        Returns
        -------
        dict
            Dictionary mapping each group label to its centroid (2d real array).
            If return_std=True, also includes key '<label>_std' with standard deviation.
        """
        z = SinusoidalEncoder.to_complex(emb)
        out = {}

        if labels is None:
            mean = np.mean(z, axis=0)
            result = np.concatenate([np.imag(mean), np.real(mean)])
            out['all'] = result
            if return_std:
                std = np.std(z, axis=0)
                out['all_std'] = np.concatenate([np.imag(std), np.real(std)])
            return out

        labels = np.asarray(labels)
        for g in np.unique(labels):
            mask = labels == g
            group_z = z[mask]
            mean = np.mean(group_z, axis=0)
            result = np.concatenate([np.imag(mean), np.real(mean)])
            out[g] = result
            if return_std:
                std = np.std(group_z, axis=0)
                out[f"{g}_std"] = np.concatenate([np.imag(std), np.real(std)])
        return out

    def verify_roundtrip(self, values, method='least_squares', scale='auto', verbose=True, return_details=False):
        """
        Perform an encode → decode → compare roundtrip and report accuracy.

        Parameters
        ----------
        values : array-like
            Original values to test.
        method : str
            Decoding method: 'least_squares', 'svd', 'optimize', or 'naive'.
        scale : float or 'auto' or None
            Scaling strategy: 'auto' uses fit_encoder(), float uses fixed scaling, None disables scaling.
        verbose : bool
            If True, prints accuracy report.
        return_details : bool
            If True, also returns the encoded array, decoded values, and residuals.

        Returns
        -------
        success : bool
            Whether all values were accurately recovered within residual tolerance.
        details : tuple, optional
            Tuple (encoded, decoded, residuals) if return_details=True
        """
        values = np.asarray(values, dtype=float).flatten()

        if scale == 'auto':
            self.fit_encoder(values)
        elif isinstance(scale, (float, int)):
            self._scale = float(scale)
            self._auto_scale_enabled = True
        elif scale is None:
            self._scale = None
            self._auto_scale_enabled = False
        else:
            raise ValueError("Invalid value for scale. Use float, 'auto', or None.")

        encoded = self.encode(values)
        decoded, residual = self.decode(encoded, method=method, return_error=True)
        decoded_array = np.asarray(decoded, dtype=float)
        error = np.abs(decoded_array - values)
        max_error = np.max(error)
        mean_error = np.mean(error)
        bad = np.where(residual > self._decode_residual_tolerance)[0]
        success = len(bad) == 0

        if verbose:
            print(f"🔁 Roundtrip verification ({method}, scale={self._scale}):")
            print(f"    Max abs error     : {max_error:.3e}")
            print(f"    Mean abs error    : {mean_error:.3e}")
            print(f"    Residual failures : {len(bad)} / {len(values)} "
                  f"(threshold = {self._decode_residual_tolerance})")
            if not success:
                print(f"    ⚠️  Failed indices: {bad.tolist()}")

        return (success, (encoded, decoded_array, residual)) if return_details else success

    @staticmethod
    def phase_unwrap(emb, normalize=False):
        """
        Perform Fourier-like phase unwrapping on sinusoidal embedding.

        Parameters
        ----------
        emb : np.ndarray
            Encoded array of shape (n, 2d).
        normalize : bool
            Whether to scale unwrapped phases to [0, 1].

        Returns
        -------
        np.ndarray
            Phase unwrapped matrix of shape (n, d), optionally normalized.
        """
        z = SinusoidalEncoder.to_complex(emb)
        phases = np.unwrap(np.angle(z), axis=0)
        if normalize:
            minv = np.min(phases, axis=0, keepdims=True)
            maxv = np.max(phases, axis=0, keepdims=True)
            return (phases - minv) / (maxv - minv + 1e-9)
        return phases


# %% Low-Level classes to store code and codefull DNA
class DNACodes(UserDict):
    """
    🧬 DNACodes
    Dictionary-like container for symbolic signal encodings at multiple scales.

    Attributes
    ----------
    meta : dict
        Metadata describing the signal and encoding parameters.
    encoded : bool
        Whether the content has been sinusoidally encoded.

    Methods
    -------
    sinencode(d_part=32, N=10000)
        Encodes symbolic segments using transformer-style sinusoidal embeddings.
    sindecode(reference_dx=None)
        Decodes sinusoidal embeddings back to symbolic segment structure.
    summary()
        Displays segment or vector counts by scale.
    plot(figsize=(12, 4), d_part=None, N=None)
        Plot method for DNACodes
    """

    def __init__(self, *args, meta=None, encoded=False, **kwargs):
        super().__init__(*args, **kwargs)
        self.meta = meta or {}
        self.encoded = encoded

    def sinencode(self, d_part=32, N=10000):
        """
        Encode symbolic segments at each scale using sinusoidal encoding grouped by letter.

        Parameters
        ----------
        d_part : int
            Number of dimensions for each component (start, width, height).
        N : int
            Frequency base for sinusoidal embedding.

        Returns
        -------
        DNACodes
            Encoded version of the current codes, grouped by letter per scale.
        """
        encoded = DNACodes(meta=self.meta.copy(), encoded=True)
        for scale, code in self.items():
            encoder = SinusoidalEncoder(d_model=3 * d_part, N=N)
            encoded[scale] = encoder.sinencode_dna_grouped(code, d_part=d_part, N=N)
        encoded.meta.update({"d_part": d_part, "N": N})
        return encoded

    def sindecode(self, reference_dx=None):
        """
        Decode sinusoidally embedded codes grouped by letter into symbolic segment structure.

        Parameters
        ----------
        reference_dx : float, optional
            Sampling interval used to reconstruct `xloc`. Defaults to meta["sampling_dt"].

        Returns
        -------
        DNACodes
            Decoded symbolic codes for each scale.
        """
        if not self.encoded:
            raise ValueError("sindecode() can only be called on encoded instances.")
        decoded = DNACodes(meta=self.meta.copy(), encoded=False)
        ref_dx = reference_dx or self.meta.get("sampling_dt", 1.0)
        ref = {"dx": ref_dx}
        d_part = self.meta.get("d_part", 32)
        N = self.meta.get("N", 10000)
        for scale, grouped in self.items():
            decoded[scale] = SinusoidalEncoder.sindecode_dna_grouped(grouped, reference_code=ref, d_part=d_part, N=N)
        return decoded

    def summary(self):
        """
        Print the number of encoded vectors or symbolic segments per scale.
        """
        for scale, content in self.items():
            if self.encoded:
                counts = {k: len(v) for k, v in content.items()}
                total = sum(counts.values())
                detail = ", ".join(f"{k}:{v}" for k, v in counts.items())
                print(f"🔢 Scale {scale} → {total} vectors ({detail})")
            else:
                print(f"🔡 Scale {scale} → {len(content['letters'])} segments")

    def plot(self, figsize=(12, 4), d_part=None, N=None):
        """
        Plot method for DNACodes: visualizes encoded vectors or symbolic segment distribution.

        Parameters
        ----------
        figsize : tuple
            Figure size for the entire plot.
        d_part : int, optional
            Number of dimensions per segment part (only for encoded).
        N : int, optional
            Frequency base (for metadata or title info).
        """
        if not self:
            print("DNACodes is empty.")
            return
        d_part = d_part or self.meta.get("d_part", 32)
        N = N or self.meta.get("N", 10000)
        if self.encoded:
            all_letters = sorted({l for scale_data in self.values() for l in scale_data})
            n_letters = len(all_letters)

            fig, axes = plt.subplots(n_letters, 3, figsize=(figsize[0], figsize[1] * n_letters), squeeze=False)
            norm = Normalize(vmin=min(self.keys()), vmax=max(self.keys()))
            colors = cm.viridis(norm(list(self.keys())))

            for i, letter in enumerate(all_letters):
                for j, field in enumerate(["start", "width", "height"]):
                    ax = axes[i][j]
                    for idx, scale in enumerate(sorted(self.keys())):
                        data = self[scale].get(letter)
                        if data is not None:
                            segment = data[:, j * d_part:(j + 1) * d_part]
                            for dim in range(min(3, d_part)):  # Plot only 3 components max
                                ax.plot(segment[:, dim], color=colors[idx], alpha=0.6)
                    if i == 0:
                        ax.set_title(field.capitalize())
                    if j == 0:
                        ax.set_ylabel(letter)
            handles = [plt.Line2D([0], [0], color=colors[i], label=f"scale {s}") for i, s in enumerate(sorted(self.keys()))]
            fig.legend(handles=handles, loc='upper right')
            fig.suptitle("Encoded DNACodes grouped by letter and field")
        else:
            # Plot bar plots for raw letter frequencies
            fig, axes = plt.subplots(1, len(self), figsize=(figsize[0], figsize[1]))
            if len(self) == 1:
                axes = [axes]
            for ax, (scale, code) in zip(axes, sorted(self.items())):
                letters = code['letters']
                unique, counts = np.unique(list(letters), return_counts=True)
                ax.bar(unique, counts)
                ax.set_title(f"Scale {scale}")
                ax.set_ylabel("Frequency")
            fig.suptitle("Raw DNACodes frequency per letter")
        plt.tight_layout()
        plt.show()
        return fig

class DNAFullCodes(dict):
    """
    🧬 DNAFullCodes(dict)
    ----------------------
    A container for symbolic full-resolution DNA-like strings or their sinusoidal embeddings,
    organized per scale.

    This structure maps each `scale` (typically corresponding to a wavelet or resolution level)
    to either:
      - a DNA-like string (`str` or `DNAstr`) representing symbolic patterns over time, or
      - a compressed embedding (dict of vectors) after sinusoidal encoding.

    It supports signal discretization, symbolic transformation, sinusoidal encoding,
    dimensionally-reduced analysis, and visual comparison of encoded motifs.

    Attributes
    ----------
    meta : dict
        Optional metadata (e.g. sampling rate, units, scale definitions, etc.).
    encoded : bool
        Whether this instance contains sinusoidally encoded data.
    unwrapped_matrix : dict, optional
        When applicable, stores a matrix {scale: ndarray (n_letters, d_model)} from compressed
        representations via `unwrap_letters_to_matrix()`.

    Methods
    -------
    sinencode(d_model=96, N=10000, operation=None)
        Encodes symbolic data with sinusoidal positional encoding. Supports per-letter reduction via
        'sum' or 'mean'. Returns a new encoded instance.

    sindecode()
        Attempts to reconstruct the symbolic string by repeating letters. Only works if the original
        operation did not compress to a single vector per letter.

    unwrap_letters_to_matrix()
        Converts compressed encodings (after 'sum' or 'mean') into (n_letters × d_model) matrices
        per scale. Required for d-space plotting.

    plot(figsize=(12, 4))
        Plots the letter-wise composition (symbolic form) or encoded means (if encoded=True).

    plot_unwrapped_matrix(figsize=(12, 4))
        Visualizes each letter’s embedding vector in d-space, with one subplot per scale.

    Example
    -------
    >>> codes = DNAFullCodes({4: 'YAABZZ'}, meta={"sampling_dt": 0.5})
    >>> encoded = codes.sinencode(operation="mean")
    >>> encoded.unwrap_letters_to_matrix()
    >>> encoded.plot_unwrapped_matrix()
    """

    def __init__(self, *args, meta=None, encoded=False, **kwargs):
        super().__init__(*args, **kwargs)
        self.meta = meta or {}
        self.encoded = encoded

    def summary(self):
        """
        Return a brief summary of the full codes per scale.

        Returns
        -------
        dict
            Mapping from scale to summary string or code length.
        """
        return {scale: f"{len(seq)} letters" for scale, seq in self.items()}

    def sinencode(self, d_model=96, N=10000, operation="sum"):
        """
        Sinencode method — encodes each letter in the DNAFullCodes as a set of sinusoidal embeddings.

        Parameters
        ----------
        d_model : int, optional
            Dimensionality of the sinusoidal embedding (default is 96).
        N : int, optional
            Maximum number of positions for the encoding (default is 10000).
        operation : str or None, optional
            If "sum", sums all encodings per letter.
            If "mean", averages all encodings per letter.
            If None, keeps the full (n_occurrences, d_model) matrix per letter.
            Raises a ValueError if the operation is not one of the above.

        Returns
        -------
        DNAFullCodes (without aggregation)
            A new DNAFullCodes instance with encoded representations and metadata.
        DNAFullCodes (with aggregation)
            A new DNAFullCodes instance with encoded and aggregated (based operation) representations and metadata.
        """
        encoded = DNAFullCodes(meta=self.meta.copy(), encoded=True)            # encoded signal
        encoded_aggregated = DNAFullCodes(meta=self.meta.copy(), encoded=True) # encoded and aggregated
        for scale, seq in self.items():
            encoder = SinusoidalEncoder(d_model=d_model, N=N)
            grouped = {ch: encoder.encode([i for i, c in enumerate(seq) if c == ch])
                       for ch in set(seq)}
            aggregated = {}
            if operation is not None:
                for ch, array in grouped.items():
                    if operation == "sum":
                        aggregated[ch] = array.sum(axis=0)
                    elif operation == "mean":
                        aggregated[ch] = array.mean(axis=0)
                    else:
                        raise ValueError(f"Unsupported operation: {operation}")
            encoded[scale] = grouped
            encoded_aggregated[scale] = aggregated
        encoded.meta.update({"d_model": d_model, "N": N, "operation": None})
        if operation:
            encoded_aggregated.meta.update({"d_model": d_model, "N": N, "operation": operation})
            encoded_aggregated.unwrap_letters_to_matrix() # we unwrap all data
        return encoded, encoded_aggregated

    def sindecode(self):
        """Sindecode method"""
        if not self.encoded:
            raise ValueError("sindecode() can only be called on encoded instances.")
        decoded = DNAFullCodes(meta=self.meta.copy(), encoded=False)
        for scale, grouped in self.items():
            letters = []
            for letter, array in grouped.items():
                letters.extend([letter] * len(array))
            decoded[scale] = DNAstr("".join(letters))
        return decoded

    def unwrap_letters_to_matrix(self):
        """
        Assemble all encoded letter vectors into a matrix of shape (n_letters, d_model)
        for each scale.

        Applies only when the encoding was performed with an operation ("sum" or "mean").
        Stores the result in self.unwrapped_matrix as a dict {scale: matrix}.
        Returns the dictionary for chaining or inspection.

        Raises
        ------
        ValueError if the encoding is not compressed (i.e., operation is None or missing),
        or if encoded entries are inconsistent in shape.
        """
        if not self.encoded:
            raise ValueError("unwrap_letters_to_matrix() can only be used on encoded instances.")

        op = self.meta.get("operation", None)
        if op not in {"sum", "mean"}:
            raise ValueError(f"Cannot unwrap: incompatible operation mode '{op}'. Must be 'sum' or 'mean'.")

        result = {}
        for scale, grouped in self.items():
            vectors = list(grouped.values())
            letters = list(grouped.keys())

            # Ensure all vectors are 1D and same length
            if not all(isinstance(v, np.ndarray) and v.ndim == 1 for v in vectors):
                raise ValueError(f"Scale {scale}: not all vectors are compressed 1D arrays.")
            d_set = {v.shape[0] for v in vectors}
            if len(d_set) > 1:
                raise ValueError(f"Scale {scale}: inconsistent vector dimensions found: {d_set}")

            mat = np.stack(vectors, axis=0)  # shape (n_letters, d)
            result[scale] = mat

        self.unwrapped_matrix = result
        return result


    def __repr__(self):
        return f"<DNAFullCodes(scales={list(self.keys())}, encoded={self.encoded})>"

    def __str__(self):
        return "\n".join([f"scale {s} → {len(v)} letters" for s, v in self.items()])

    def plot(self, figsize=(12, 4)):
        """
        Plot method for DNAFullCodes: visualizes encoded vectors or DNA string composition.

        Parameters
        ----------
        figsize : tuple
            Figure size for the entire plot.
        """
        if not self:
            print("DNAFullCodes is empty.")
            return
        if self.encoded:
            all_letters = sorted({ch for self in self.values() for ch in self})
            n_letters = len(all_letters)
            scales = sorted(self.keys())
            fig, axs = plt.subplots(n_letters, 1, figsize=(figsize[0], figsize[1]*n_letters), squeeze=False)
            axs = axs[:, 0]
            cmap = cm.get_cmap('viridis', len(scales))
            for i, letter in enumerate(all_letters):
                ax = axs[i]
                for j, scale in enumerate(scales):
                    scale_data = self[scale]
                    if letter not in scale_data:
                        continue
                    arr = scale_data[letter]
                    x = np.arange(arr.shape[1])
                    y_mean = arr.mean(axis=0)
                    ax.plot(x, y_mean, label=f"scale {scale}", color=cmap(j))
                ax.set_title(f"Letter '{letter}'")
                ax.set_ylabel("Mean Encoding")
                ax.set_xlabel("Embedding Dimension")
                ax.grid(True)
                ax.legend()

            plt.tight_layout()
            plt.show()
            fig.suptitle("Encoded DNAFullCodes grouped by letter")
        else:
            all_letters = sorted({c for v in self.values() if not isinstance(v, dict) for c in str(v)})
            n_letters = len(all_letters)
            fig, axes = plt.subplots(1, len(self), figsize=(figsize[0], figsize[1]))
            if len(self) == 1:
                axes = [axes]
            for ax, (scale, seq) in zip(axes, sorted(self.items())):
                unique, counts = np.unique(list(str(seq)), return_counts=True)
                ax.bar(unique, counts)
                ax.set_title(f"Scale {scale}")
                ax.set_ylabel("Frequency")
            fig.suptitle("Raw DNAFullCodes frequency per letter")
        plt.tight_layout()
        plt.show()
        return fig

    def plot_unwrapped_matrix(self, figsize=(12, 4)):
        """
        Plot each letter's encoded vector in the abstract embedding space (d-space).
        One curve per letter, one subplot per scale.

        Requires `unwrap_letters_to_matrix()` to have been called.

        Parameters
        ----------
        figsize : tuple
            Base figure size. Height will be scaled based on number of scales.

        Returns
        -------
        matplotlib.figure.Figure
            The generated matplotlib figure.
        """
        if not hasattr(self, "unwrapped_matrix"):
            raise AttributeError("unwrapped_matrix not found. Run unwrap_letters_to_matrix() first.")
        scales = sorted(self.unwrapped_matrix.keys())
        n_scales = len(scales)
        fig, axs = plt.subplots(n_scales, 1, figsize=(figsize[0], figsize[1]*n_scales), squeeze=False)
        axs = axs[:, 0]

        for i, scale in enumerate(scales):
            matrix = self.unwrapped_matrix[scale]
            letters = list(self[scale].keys())
            for j, vec in enumerate(matrix):
                axs[i].plot(vec, label=f"'{letters[j]}'")
            axs[i].set_title(f"Scale {scale}")
            axs[i].set_xlabel("Encoding Dimension")
            axs[i].set_ylabel("Value")
            axs[i].legend()
            axs[i].grid(True)

        fig.suptitle("Encoded letter vectors per scale in d-space")
        plt.tight_layout()
        return fig


# %% Main classes DNAsignal, DNAstr, DNApairwiseAnalysis

# ------------------------
# DNAsignal
# ------------------------
class DNAsignal:
    """
    DNAsignal(signal):
    ==================
    A class to encode a numerical signal (typically a 1D GC-MS trace, NMR/FTIR/Raman spectra
    or time series) into a DNA-like symbolic representation using wavelet analysis.
    This symbolic coding enables fast comparison, search, and alignment of signal features
    using abstracted patterns (e.g., 'YAZB').

    The class supports:
    - Continuous Wavelet Transform (CWT) using Ricker wavelets
    - Symbolic conversion of wavelet features to DNA-like letters (A, B, C, X, Y, Z)
    - Visualization of CWTs, symbolic encodings, and signal overlays
    - Reversible decoding of symbolic segments back into approximate signals
    - Substring extraction and matching
    - Storage of multi-scale representations (multi-resolution DNA encoding)

    This class is part of the symbolic signal transformation pipeline (`sig2dna`), compatible with `signal`, `peaks`, and `signal_collection`.

    Parameters
    ----------
    signal : signal
        An instance of the `signal` class representing a sampled waveform. It must
        have `x`, `y`, and a valid name or identifier.
    encode : bool (default = False)
        Launch encoders = ["compute_cwt","encode_dna","encode_dna_full"] if True
    plot : bool (default = False)
        Plot with plotter = ["plot_signals","plot_transforms","plot_codes"]

    Attributes
    ----------
    signal : np.ndarray
        Original numerical signal (values only; x stored via `sampling_dt`).
    dtype : data-type
        Data type of the stored signal.
    sampling_dt : float
        Sampling interval along the x-axis.
    dx : float (depreciated)
        The nominal x-resolution of the signal (automatically derived).
    n : int
        Number of points in the signal (length of x/y arrays).
    name : str
        Name of the signal.
    x_label : str
        Label of the x-axis.
    x_unit : str
        Unit of the x-axis.
    y_label : str
        Label of the y-axis.
    y_unit : str
        Unit of the y-axis.
    scales : list[int]
        List of scales used for encoding.
    codes : DNACodes
        Dictionary-like container of symbolic triplet encodings by scale.
    codesfull : DNAFullCodes
        Dictionary-like container of full-resolution symbolic strings per scale.
    scales : array-like
        Set of scales used in the Continuous Wavelet Transform (powers of 2 by default).
    transforms : signal_collection
        Stores the CWT-transformed signals, each as a `signal` object, indexed by scale.
    codes : dict[int, DNAstr]
        Symbolic codes by scale level. Each is a `DNAstr` object representing
        the symbolic sequence at that scale.
    codesfull : dict[int, DNAstr]
        Same as `codes`, but uses full resolution symbolic representation.
    sincodesfull : dict[int, DNAFullCodes]
        Sinusoidal position encoded DNAstr (without aggregation)
    sincodesfull_aggregated
        Sinusoidal position encoded DNAstr (with aggregation)
    peaks : peaks
        Optional `peaks` object used to index real peak positions from the signal.
    codebook : dict
        Mapping between symbolic characters (A, B, C, X, Y, Z, _) and wavelet features.
    generator : str
        Name of the wavelet basis used (default: 'ricker').

    Methods
    -------
    normalize_signal(mode="zscore+shift")
        Normalizes the internal signal (preserves positivity).
    compute_cwt(scales=None, normalize=False)
        Computes the Continuous Wavelet Transform using the Ricker wavelet
        and stores the transformed signals in `transforms`.
    sparsify_cwt(self, scale: Union[int, float], threshold: float, inplace: bool = True)
        Zero out wavelet coefficients below a threshold for a specific scale.
    encode_dna()
        Encodes each scale's transformed signal into a symbolic `DNAstr` sequence
        using local maximum coding (ABCXYZ).
    encode_dna_full()
        Encodes signals at each scale using the full encoding scheme, preserving
        the flat regions (_) and finer symbolic transitions.
    plot_signals()
        Plots signals.
    plot_codes(scale)
        Plots both the wavelet transform and the symbolic code for the given scale.
    plot_transforms
        Plots the stored CWT-transformed signals as a signal collection.
    plot_scalogram():
        Plots a scalogram with two subplots
    decode_dna(scale)
        Reconstructs the approximate signal for a given scale from its DNA encoding.
    __getitem__(scale)
        Shortcut to access the `DNAstr` object for a specific scale.
    summary()
        Returns a dictionary summarizing encoded scales and metadata.
    has(scale)
        Checks whether a DNAstr encoding exists at the given scale.
    pseudoinverse(scales=None, rank=None, return_weights=False, name=None):
        Approximates signal reconstruction via pseudo-inverse using stored CWT coefficients

    Static pairwise distance methods
    --------------------------------
    _pairwiseEntropyDistance(list of DNAstr objects, scale)
        Return a DNApairwiseAnalysis instance based on the mutually exclusive information
        after DNA/code alignment.
    _pairwiseJaccardMotifDistance(list of DNAstr objects, scale)
        Return a DNApairwisedistance based on the presence/absence of a pattern (default=YAZB)
    _pairwiseJensenShannonDistance(list of DNAstr objects, scale)
        Return a DNApairwisedistance based on the Jensen-Shannon distance of a pattern
    _pairwiseLevenshteinDistance(list of DNAstr objects, scale)
        Return a DNApairwisedistance based on the Levenshtein Distance

    Methods: Symbolic Encoding
    --------------------------
    encode_dna(scales=None)
        Convert signal into triplet-based symbolic encoding (per scale).
    encode_dna_full(scales=None, resolution='index', repeat=True, n_points=None)
        Generate full-resolution DNA strings by repeating letters (to `codesfull`).

    Methods: Sinusoidal Encoding
    ----------------------------
    sinencode_dna(scales=None, d_part=32, N=10000)
        Encode symbolic segments (from `codes`) into sinusoidal vectors.
    sinencode_dna_full(scales=None, d_part=32, N=10000)
        Encode symbolic full strings (from `codesfull`) into sinusoidal vectors.

    Methods: Sinusoidal Decoding (Static)
    -------------------------------------
    sindecode_dna(grouped_embeddings, reference_dx=1.0, d_part=32, N=10000)
        [static] Decode grouped sinusoidal embeddings → `DNACodes` structure.
    sindecode_dna_full(grouped_embeddings, reference_dx=1.0, d_part=32, N=10000)
        [static] Decode full sinusoidal embeddings → `DNAFullCodes` structure.

    Methods: Signal Reconstruction
    ------------------------------
    tosignal(scale=None, codes_attr='codes')
        Reconstruct approximate signal (as `signal` instance) from symbolic encoding.

    Examples
    --------
    >>> S = signal.from_peaks(...)  # define a signal
    >>> dna = DNAsignal(S)
    >>> dna.compute_cwt()
    >>> dna.encode_dna()
    >>> dna.codes[4]
    DNAstr("AAAZZZYY...")
    >>> dna.plot_codes(4)
    >>> dna.codes[4].find("YAZB")
    """

    def __init__(self, signal_obj, sampling_dt=1.0, dtype=np.float64,
                 encode=False,
                 encoder=["compute_cwt","encode_dna","encode_dna_full"],
                 scales=[1,2,4,8,16,32],
                 x_label="index", x_unit = "-", y_label="Intensity", y_unit = "",
                 plot=False,
                 plotter=["plot_signals","plot_transforms","plot_codes"]):
        """
        Initialize DNAsignal with a signal object or 1D array.

        Parameters
        ----------
        signal_obj : signal or np.ndarray
            The input signal (preferably a signal object).
        sampling_dt : float
            Sampling interval (used for filtering).
        dtype : data-type
            Data type to cast the signal (default: np.float64).
        """
        if isinstance(signal_obj, signal):
            self.signal_obj = signal_obj.copy()
            self.signal = np.asarray(signal_obj.y, dtype=dtype)
            self.x = signal_obj.x
            self.x_label = signal_obj.x_label or x_label
            self.x_unit = signal_obj.x_unit or x_unit
            self.y_label = signal_obj.y_label or y_label
            self.y_unit = signal_obj.y_unit or y_unit
            self.name = signal_obj.name
            self.sampling_dt = signal_obj.x[1] - signal_obj.x[0]
        else:
            self.signal_obj = None
            self.signal = np.asarray(signal_obj, dtype=dtype)
            self.x = np.arange(len(self.signal)) * sampling_dt
            self.x_label = x_label
            self.x_unit =  x_unit
            self.y_label = y_label
            self.y_unit = y_unit
            self.name = "array"
            self.sampling_dt = sampling_dt
        self.dtype = dtype
        self.filtered_signal = None
        self.scales = []
        self.cwt_coeffs = {}  # scale -> array
        self.codes = DNACodes()  # dict-like[scale]-> {letters, widths, heights}

        # encode and plots on request
        if plot and "plot_signals" in plotter:
            self.plot_signals()
        if encode:
            self.compute_cwt(scales=scales) # minimum encoder
            if plot and "plot_transforms" in plotter:
                self.plot_transforms()
            if ("encode_dna" in encoder) or ("encode_dna_full" in encoder):
                self.encode_dna(scales=scales)
                if plot and "plot_codes" in plotter:
                    self.plot_codes(scales=scales)
            if ("encode_dna_full" in encoder):
                self.encode_dna_full()

    def normalize_signal(self, mode="zscore+shift"):
        """
        Normalize the internal signal using one of several strategies that ensure positivity.

        Parameters
        ----------
        mode : str
            Normalization mode passed to `signal.normalize()`.
            See `signal.normalize()` for available modes.

        Raises
        ------
        AttributeError
            If `signal` attribute is missing or of the wrong type.
        """
        if not hasattr(self, 'signal') or not isinstance(self.signal, signal):
            raise AttributeError("DNAsignal must have a valid `signal` attribute of type `signal`.")
        self.signal.normalize(mode=mode, inplace=True)

    def has(self, scale):
        """
        Check if a DNA encoding exists for the specified scale.

        Parameters
        ----------
        scale : int
            The wavelet scale to check.

        Returns
        -------
        bool
            True if a symbolic DNAstr encoding exists at the given scale, False otherwise.

        Examples
        --------
        >>> dna.has(4)
        True
        >>> dna.has(16)
        False
        """
        return scale in self.codes

    @staticmethod
    def apply_baseline_filter(signal, w=None, k=2, delta_t=1.0):
        """
        Apply baseline filtering using moving median and local Poisson-based thresholding.

        Parameters
        ----------
        signal : np.ndarray
            Input signal (expected to be non-negative or baseline-dominated).
        w : int or None
            Window size for baseline and statistics (must be odd).
            Defaults to max(11, 2% of signal length).
        k : float
            Bienaymé-Tchebychev multiplier.
        delta_t : float
            Sampling time step.

        Returns
        -------
        filtered : np.ndarray
            Signal with baseline removed and low-intensity noise suppressed.

        Note
        ----
        This method is static, use signal.apply_baseline_filter() whenever appropriate instead.

        """
        signal = np.asarray(signal)
        m = len(signal)

        # Determine appropriate window width
        if w is None:
            w = max(11, int(0.01 * m))
        if w % 2 == 0:
            w += 1
        if w >= m:
            raise ValueError(f"Window width w={w} must be smaller than signal length {m}.")

        # Step 1: remove baseline via moving median
        baseline = medfilt(signal, kernel_size=w)
        s = signal - baseline
        s[s < 0] = 0  # force non-negativity

        # Step 2: moving mean and std (uniform filter = moving average)
        # signal.apply_baseline_filter() uses np.sliding_window_view() instead.
        mean = uniform_filter1d(s, size=w, mode='nearest')
        sq = uniform_filter1d(s**2, size=w, mode='nearest')
        std = np.sqrt(np.maximum(sq - mean**2, 0))

        # Step 3: Poisson λ from cv = std / mean
        with np.errstate(divide='ignore', invalid='ignore'):
            cv = np.where(mean > 0, std / mean, 0)
            lam = np.where(cv > 0, 1 / (cv**2), 0)

        # Step 4: BT thresholding
        threshold = k * np.sqrt(10 * lam * delta_t)
        s[s < threshold] = 0
        return s

    def compute_cwt(self, scales=None, apply_filter=False, wavelet="mexh"):
        """
        Compute Continuous Wavelet Transform (CWT) using the Mexican Hat wavelet.

        Parameters
        ----------
        scales : list, int, or None
            List of scales (or a single scale) to compute. If None, default to [1, 2, 4, 8, 16].
        apply_filter : bool
            Whether to apply a baseline filter to the input signal before transforming.

        wavelet : str (default='mexh')
            The name of the PyWavelets-compatible wavelet.
        Sets
        ----
        self.scales : list
            The list of actual scales used.
        self.filtered_signal : ndarray
            Filtered or raw signal used for CWT.
        self.cwt_coeffs : dict
            Dictionary mapping each scale to its 1D coefficient array.
        self.transforms : signal_collection
            Collection of `signal` objects storing the transformed signals for each scale.
        """
        if scales is None:
            scales = [2 ** i for i in range(5)]
        if not isinstance(scales,(list,tuple)):
            scales = [scales]
        self.scales = scales
        if apply_filter:
            self.filtered_signal = self.apply_baseline_filter(self.signal,delta_t=self.sampling_dt)
        else:
            self.filtered_signal = self.signal
        self.transforms = signal_collection()
        for scale in self.scales:
            coef, _ = pywt.cwt(self.filtered_signal, [scale], wavelet)
            self.cwt_coeffs[scale] = coef[0]  # Access directly via scale
            sig = signal(
                x=self.x.copy(),
                y=coef[0],
                name=f"CWT_scale_{scale}",
                x_label="x",
                y_label="CWT amplitude",
                y_unit="a.u.",
                source="CWT")
            self.transforms.append(sig)

    def pseudoinverse(self, scales=None, rank=None, return_weights=False, name=None):
        """
        Approximate signal reconstruction via pseudo-inverse using stored CWT coefficients.

        Parameters
        ----------
        scales : list, float, int, or None
            Scales to include in the reconstruction. If None, all scales in self.cwt_coeffs are used.
        rank : int or None
            Optional truncation rank for the SVD decomposition (for denoising or dimensionality reduction).
        return_weights : bool
            If True, also return the weights (contributions) of each scale.
        name : str or None
            Optional name for the returned signal. Defaults to "pseudoinverse" with included scales.

        Returns
        -------
        reconstructed_signal : signal
            Reconstructed signal instance from the pseudo-inverse of the CWT decomposition.
        weights : np.ndarray, optional
            Returned only if `return_weights=True`, gives the contribution of each scale.

        Raises
        ------
        ValueError
            If CWT coefficients are not available.
        """
        if not hasattr(self, "cwt_coeffs") or not self.cwt_coeffs:
            raise ValueError("No CWT coefficients found. Run compute_cwt() first.")

        if scales is None:
            scales = sorted(self.cwt_coeffs.keys())
        elif isinstance(scales, (int, float)):
            scales = [scales]

        # Assemble wavelet coefficient matrix
        W = np.vstack([self.cwt_coeffs[s] for s in scales])  # shape: (n_scales, n_points)

        # SVD-based pseudoinverse
        U, S, Vt = np.linalg.svd(W, full_matrices=False)
        if rank is not None:
            U, S, Vt = U[:rank], S[:rank], Vt[:rank]

        recon_proj = np.dot(U.T, W)
        weights = np.dot(np.diag(1 / S), recon_proj)
        recon_y = np.dot(weights.T, U).sum(axis=1)

        # Build signal object
        sig_name = name or f"pseudoinverse_{','.join(map(str, scales))}"
        recon_signal = signal(
            x=self.x.copy(),
            y=recon_y,
            name=sig_name,
            type="reconstructed",
            x_label=self.signal.x_label,
            x_unit=self.signal.x_unit,
            y_label=self.signal.y_label,
            y_unit=self.signal.y_unit,
            source=f"pseudo-inverse from scales {scales}",
            color="black",
            linestyle="--",
            message=f"Reconstructed using pseudo-inverse with rank={rank}"
        )
        return (recon_signal, weights) if return_weights else recon_signal


    def sparsify_cwt(self, scale=None, threshold=None, inplace=True):
        """
        Sparsify CWT coefficients by zeroing values below a threshold.

        Parameters
        ----------
        scale : int, float, list, or None
            Scale(s) to sparsify. If None, all available scales in self.cwt_coeffs are used.
        threshold : float or None
            Absolute value below which coefficients are set to zero.
            If None, uses 1% of the maximum absolute value at each scale.
        inplace : bool
            If True, modifies current instance. If False, returns a modified copy.

        Returns
        -------
        DNAsignal or None
            Modified copy if inplace is False, otherwise None.

        Raises
        ------
        ValueError
            If scale(s) not found in self.cwt_coeffs.
        """
        if not hasattr(self, "cwt_coeffs") or not self.cwt_coeffs:
            raise ValueError("No CWT coefficients found. Run compute_cwt() first.")

        # Normalize scale input
        if scale is None:
            target_scales = list(self.cwt_coeffs.keys())
        elif isinstance(scale, (int, float)):
            target_scales = [scale]
        elif isinstance(scale, (list, tuple)):
            target_scales = scale
        else:
            raise TypeError(f"Invalid type for scale: {type(scale)}")

        # Check all requested scales exist
        missing = [s for s in target_scales if s not in self.cwt_coeffs]
        if missing:
            raise ValueError(f"Scales not found in computed CWT coefficients: {missing}")

        # Work in-place or copy
        target = self if inplace else self.copy()

        for s in target_scales:
            coeffs = target.cwt_coeffs[s]
            thres = threshold if threshold is not None else 0.01 * np.max(np.abs(coeffs))
            target.cwt_coeffs[s] = np.where(np.abs(coeffs) < thres, 0.0, coeffs)

            if hasattr(target, "transforms"):
                for sig in target.transforms:
                    if sig.name == f"CWT_scale_{s}":
                        sig.y = target.cwt_coeffs[s]
                        sig._events("sparsify", {"scale": s, "threshold": thres})
                        break

        return None if inplace else target


    # --------------------------[ DNA coding segment-wise (no repetitions) ]----------------------------

    def encode_dna(self, scales=None):
        """
        Encode each transformed signal into a symbolic DNA-like sequence of monotonic segments.

        Parameters
        ----------
        scales : list, int, or None
            List of scales (or a single scale) to encode. If None, use self.scales.

        The encoding detects strictly monotonic (or flat) segments and labels them with symbolic letters:
        - A: crosses 0 upward (neg → pos)
        - Z: crosses 0 downward (pos → neg)
        - B: strictly increasing negative segment
        - Y: strictly decreasing negative segment
        - C: strictly increasing positive segment
        - X: strictly decreasing positive segment
        - _: flat or ambiguous segment

        Sets
        ----
        self.codes : dict
            Dictionary mapping each scale to a struct with:
                - letters : str (symbolic encoding)
                - widths  : list of float (x-span of each segment)
                - heights : list of float (y-delta of each segment)
                - iloc    : list of index-pair tuples (start, end+1)
                - xloc    : list of x-span tuples (x_start, x_end)
                - dx      : segment step (dx)
        """
        if scales is None:
            scales = self.scales
        if not isinstance(scales,(list,tuple)):
            scales = [scales]
        if not hasattr(self, 'cwt_coeffs') or not self.cwt_coeffs:
            self.compute_cwt(scales)
        dx = self.x[1]-self.x[0]
        for scale in scales:
            coef = self.cwt_coeffs[scale]
            letters, widths, heights, iloc, xloc = [], [], [], [], []
            monotonic = np.diff(coef)
            mono_sign = np.sign(monotonic)
            sign_changes = np.where(np.diff(mono_sign) != 0)[0] + 1  # +1 because diff shortens by 1
            start_idx = 0
            segment_ends = np.append(sign_changes, len(coef) - 1)
            for count, idx in enumerate(segment_ends):
                xsegment = self.x[start_idx:idx + 1]
                segment = coef[start_idx:idx + 1]
                if len(segment) < 2:
                    letters.append('_')
                else:
                    start, end = segment[0], segment[-1]
                    letter = self._get_letter(start, end)
                    letters.append(letter)
                widths.append(xsegment[-1] - xsegment[0])
                heights.append(segment[-1] - segment[0])
                xloc.append((xsegment[0],xsegment[-1]))
                is_last = count == len(segment_ends) - 1
                iloc.append((start_idx, idx + 1 if is_last else idx))
                start_idx = idx #idx + 1 (segments are continguous)
            # update codes
            self.codes[scale] = {'letters': ''.join(letters),
                                 'widths': widths,
                                 'heights': heights,
                                 'xloc':xloc,
                                 'iloc':iloc,
                                 'dx':dx}
            return self.codes # for chaining

    def sinencode_dna(self, scales=None, d_part=32, N=10000):
        """
        Encode `self.codes` into sinusoidal embeddings (grouped by letter).

        Sets:
            - self.codes (DNACodes): Encoded version of original codes.
        """
        if not isinstance(self.codes, DNACodes):
            raise ValueError("self.codes must be a DNACodes instance")
        if scales is None:
            scales = self.scales
        codes = DNACodes({s: self.codes[s] for s in scales if s in self.codes},
                                  meta=self.codes.meta, encoded=False)
        self.sincodes = codes.sinencode(d_part=d_part, N=N)
        return self.sincodes # for chaining

    @staticmethod
    def sindecode_dna(grouped_embeddings, reference_dx=1.0, d_part=32, N=10000):
        """
        Decode sinusoidal grouped embeddings into a DNACodes structure.

        Parameters
        ----------
        grouped_embeddings : dict
            {scale: {letter: np.ndarray}} sinusoidal representations
        reference_dx : float
            Sampling resolution used to reconstruct xloc and iloc
        d_part : int
            Dimensionality per component (start, width, height)
        N : int
            Frequency base

        Returns
        -------
        DNACodes
            Decoded symbolic structure
        """
        decoded = DNACodes(meta={"sampling_dt": reference_dx}, encoded=False)
        for scale, grouped in grouped_embeddings.items():
            decoded[scale] = SinusoidalEncoder.sindecode_dna_grouped(
                grouped, reference_code={"dx": reference_dx}, d_part=d_part, N=N
            )
        return decoded


    # --------------------------[ coding in full letters (with repetitions) ]----------------------------

    def encode_dna_full(self, scales=None, resolution='index', repeat=True, n_points=None):
        """
        Convert symbolic codes into DNA-like strings by repeating letters proportionally to their span.

        Parameters
        ----------
        scales : list, int, or None
            List of scales (or a single scale) to convert. If None, use self.scales.
        resolution : {'index', 'x'}
            Repetition mode:
                - 'index': repeat letters by number of indices (j - i from iloc)
                - 'x'    : interpolate letter values over physical x-axis distance (xloc)
        repeat : bool
            If True, repeat or interpolate letters to form a string of desired resolution.
            If False, return the symbolic sequence without repetition.
        n_points : int or None
            Used only for resolution='x' to control the number of interpolation points.
            If None, defaults to ~10 points per x-unit.

        Returns
        -------
        dict
            Dictionary mapping each scale to its DNA-like string.

        Sets
        ----
        self.codesfull : dict
            Dictionary storing the resulting full DNA-like string per scale.
        """
        if scales is None:
            scales = self.scales
        elif isinstance(scales, (int, float)):
            scales = [scales]
        elif not isinstance(scales, (list, tuple)):
            raise TypeError("scales must be a list, int, or None")

        if not hasattr(self, 'codes') or not self.codes:
            self.encode_dna(scales)

        result = DNAFullCodes()
        for scale in scales:
            if scale not in self.codes:
                self.encode_dna([scale])
            code = self.codes[scale]
            if not repeat:
                result[scale] = code['letters']
                continue
            if resolution == 'index':
                sequence = ''.join(
                    letter * max(1, (j - i)) for letter, (i, j) in zip(code['letters'], code['iloc'])
                )
            elif resolution == 'x':
                if n_points is None:
                    total_span = sum(x2 - x1 for x1, x2 in code['xloc'])
                    n_points = int(np.ceil(total_span * 10))  # adjustable density
                x_start = code['xloc'][0][0]
                x_end = code['xloc'][-1][1]
                grid = np.linspace(x_start, x_end, n_points, endpoint=True)
                centers = [(x1 + x2) / 2 for x1, x2 in code['xloc']]
                # Encode letters as integers
                unique_letters = list(OrderedDict.fromkeys(code['letters']))
                letter_to_int = {ch: i for i, ch in enumerate(unique_letters)}
                int_to_letter = {i: ch for ch, i in letter_to_int.items()}

                int_vals = [letter_to_int[ch] for ch in code['letters']]
                interp_func = interp1d(centers, int_vals, kind='nearest',
                                       bounds_error=False, fill_value=int_vals[-1])
                interpolated = np.round(interp_func(grid)).astype(int)
                sequence = ''.join(int_to_letter[i] for i in interpolated)
            else:
                raise ValueError("resolution must be either 'index' or 'x'")
            result[scale] = DNAstr(sequence,
                                   dx=code["dx"],
                                   iloc=(code["iloc"][0][0],code["iloc"][-1][-1]),
                                   xloc=(code["xloc"][0][0],code["xloc"][-1][-1]),
                                   x_label=self.x_label, x_unit=self.x_unit)
        self.codesfull = result
        return self.codesfull # for chaining

    def sinencode_dna_full(self, d_model=96, N=10000, operation=None):
        """
        🌀 Encode full-resolution DNA-like strings into sinusoidal embeddings grouped by letter.

        Parameters
            ----------
        d_model : int, optional
            Dimensionality of the sinusoidal embedding (default is 96).
        N : int, optional
            Maximum number of positions for the encoding (default is 10000).
        operation : str or None, optional
            If "sum", sums all encodings per letter.
            If "mean", averages all encodings per letter.
            If None, keeps the full (n_occurrences, d_model) matrix per letter.
            Raises a ValueError if the operation is not one of the above.

        Sets
        ----
        self.codesfull : DNAFullCodes (if not already set)
            Full-resolution symbolic strings at each scale.
        self.codesfull_encoded : DNAFullCodes
            Sinusoidally encoded version of the full DNA strings.

        """
        if not hasattr(self, 'codesfull') or not self.codesfull:
            self.encode_dna_full()  # generate codesfull from codes
        if not isinstance(self.codesfull, DNAFullCodes):
            raise TypeError("self.codesfull must be an instance of DNAFullCodes.")
        # Attach metadata if missing
        self.codesfull.meta.update({
            "sampling_dt": self.sampling_dt,
            "d_model": d_model,
            "N": N
        })
        self.sincodesfull,self.sincodesfull_aggregated = \
            self.codesfull.sinencode(d_model=d_model, N=N, operation=operation)
        return self.sincodesfull # for chaining


    @staticmethod
    def sindecode_dna(grouped_embeddings, reference_dx=1.0, d_part=32, N=10000):
        """
        Decode sinusoidal grouped embeddings into a DNACodes structure.

        Parameters
        ----------
        grouped_embeddings : dict
            {scale: {letter: np.ndarray}} sinusoidal representations
        reference_dx : float
            Sampling resolution used to reconstruct xloc and iloc
        d_part : int
            Dimensionality per component (start, width, height)
        N : int
            Frequency base

        Returns
        -------
        DNACodes
            Decoded symbolic structure
        """
        decoded = DNACodes(meta={"sampling_dt": reference_dx}, encoded=False)
        for scale, grouped in grouped_embeddings.items():
            decoded[scale] = SinusoidalEncoder.sindecode_dna_grouped(
                grouped, reference_code={"dx": reference_dx}, d_part=d_part, N=N
            )
        return decoded



    # --------------------------[ back to signal ]----------------------------
    def tosignal(self, scale=None, codes_attr="codes"):
        """
        Reconstruct an approximate signal from symbolic encodings.

        Parameters
        ----------
        scale : int or None
            Scale level to use (defaults to first available if None).
        codes_attr : str
            Attribute from which to decode ('codes' or 'codesfull').

        Returns
        -------
        signal
            An approximate `signal` object reconstructed from symbolic information.
        """
        from sig2dna3.signal import signal  # Adjust path to your actual signal module

        if not hasattr(self, codes_attr):
            raise AttributeError(f"DNAsignal has no attribute '{codes_attr}'")

        codes = getattr(self, codes_attr)
        if not codes:
            raise ValueError(f"No symbolic data found in '{codes_attr}'")

        if scale is None:
            scale = next(iter(codes))

        code = codes[scale]
        if isinstance(code, str):  # from codesfull
            raise TypeError("Cannot reconstruct signal from 'codesfull' directly. Use decode + pattern extraction.")

        x_centers = [(x1 + x2) / 2 for x1, x2 in code["xloc"]]
        heights = code["heights"]

        return signal(
            x=np.array(x_centers),
            y=np.array(heights),
            name=f"decoded@{scale}",
            type="synthetic",
            x_label=self.x_label,
            x_unit=self.x_unit,
            y_label=self.y_label,
            y_unit=self.y_unit,
            source="decoded"
        )


    # --------------------------[ low-level methods ]----------------------------
    @staticmethod
    def _get_letter(start, end, tol=1e-12):
        """
        Determine letter based on monotonicity and signal range.

        Returns one of:
            - 'A': Negative to Positive (crosses zero upward)
            - 'Z': Positive to Negative (crosses zero downward)
            - 'B': Increasing Negative (concave up)
            - 'Y': Decreasing Negative (concave down)
            - 'C': Increasing Positive (concave up)
            - 'X': Decreasing Positive (concave down)
            - '_': Flat or undefined
        """
        # Optional tolerance for numerical 0
        s = start if abs(start) > tol else 0.0
        e = end   if abs(end)   > tol else 0.0
        if s == e:
            return '_'
        # Zero-crossings
        if s < 0 and e > 0:
            return 'A'
        if s > 0 and e < 0:
            return 'Z'
        # Entirely negative (possibly with 0 at one end)
        if s <= 0 and e <= 0:
            return 'B' if e > s else 'Y'
        # Entirely positive (possibly with 0 at one end)
        if s >= 0 and e >= 0:
            return 'C' if e > s else 'X'
        return '_'

    @staticmethod
    def _get_triangle_from_letter(letter,x0,y0,w,h):
        """returns the triangle (counter-clockwise, x are incr) """
        if letter == "A":
            x = (x0,x0+w,x0+w)
            y = (y0,y0,y0+h)
        elif letter == "B":
            x = (x0,x0+w,x0)
            y = (y0,y0+h,y0+h)
        elif letter == "C":
            x = (x0,x0+w,x0+w)
            y = (y0,y0,y0+h)
        elif letter =="X":
            x = (x0,x0,x0+w)
            y = (y0,y0+h,y0+h)
        elif letter == "Y":
            x = (x0,x0+w,x0+w)
            y = (y0,y0,y0+h)
        elif letter == "Z":
            x = (x0,x0,x0+w)
            y = (y0,y0+h,y0+h)
        elif letter == "_":
            x = (x0,x0+w,x0+w)
            y = (y0,y0,y0)
        else:
            raise ValueError(f"the letter {letter} is unknown/undocumented")
        return x,y


    @property
    def _is_letter_crossing(self):
        """Return True if the letter codes for a segment crossing y=0"""
        return {"A":True,"B":False,"C":False,"X":False,"Y":False,"Z":True,"_":False}
    @property
    def _is_letter_crossing_from_positive(self):
        """Return True if the letter codes for a segment crossing y=0 from y>0"""
        return {"A":False,"B":False,"C":False,"X":False,"Y":False,"Z":True,"_":False}
    @property
    def _is_letter_crossing_from_negative(self):
        """Return True if the letter codes for a segment crossing y=0 from y<0"""
        return {"A":True,"B":False,"C":False,"X":False,"Y":False,"Z":False,"_":False}
    @property
    def _is_letter_increasing(self):
        """Return True if the letter codes for an increasing segment"""
        return {"A":True,"B":True,"C":True,"X":False,"Y":False,"Z":False,"_":False}
    @property
    def _is_letter_decreasing(self):
        """Return True if the letter codes for an decreasing segment"""
        return {"A":False,"B":False,"C":False,"X":True,"Y":True,"Z":True,"_":False}
    @property
    def _is_letter_constant(self):
        """Return True if the letter codes for a constant segment, y=0"""
        return {"A":False,"B":False,"C":False,"X":False,"Y":False,"Z":False,"_":True}
    @property
    def _is_letter_starting_positive(self):
        """Return True if the letter codes for a segment starting with y>0"""
        return {"A":False,"B":False,"C":True,"X":True,"Y":False,"Z":True,"_":False}
    @property
    def _is_letter_starting_negative(self):
        """Return True if the letter codes for a segment crossing y=0"""
        return {"A":True,"B":True,"C":False,"X":False,"Y":True,"Z":False,"_":False}
    @property
    def _is_letter_ending_positive(self):
        """Return True if the letter codes for a segment ending with y>0"""
        return {"A":True,"B":False,"C":True,"X":True,"Y":False,"Z":False,"_":False}
    @property
    def _is_letter_ending_negative(self):
        """Return True if the letter codes for a segment ending with y<0"""
        return {"A":False,"B":True,"C":False,"X":False,"Y":True,"Z":True,"_":False}
    @property
    def _is_letter_starting_from_zero(self):
        """Return True if the letter codes for a segment staring from 0"""
        return {"A":False,"B":False,"C":True,"X":False,"Y":True,"Z":False,"_":False}
    @property
    def _is_letter_ending_at_zero(self):
        """Return True if the letter codes for a segment ending at 0"""
        return {"A":False,"B":True,"C":False,"X":True,"Y":False,"Z":False,"_":False}


    def get_entropy(self, scale):
        """Calculate Shannon entropy for encoded signal."""
        letters = self.codes[scale]['letters']
        _, counts = np.unique(list(letters), return_counts=True)
        return entropy(counts, base=2)

    def get_code(self, scale):
        """Retrieve encoded data for a specific scale."""
        return self.codes[scale]

    def find_sequence(self, pattern, scale):
        """Find occurrences of a specific letter pattern in encoded sequence."""
        sequence = self.codes[scale]['letters']
        positions = []
        idx = sequence.find(pattern)
        while idx != -1:
            positions.append(idx)
            idx = sequence.find(pattern, idx + 1)
        return positions

    def reconstruct_signal(self, scale, return_signal=True):
        """
        Reconstruct the signal from symbolic features (e.g., YAZB).

        Parameters
        ----------
        scale : int
            Scale to use for reconstruction.
        return_signal : bool
            If True, return a `signal` object. Else return y array.

        Returns
        -------
        signal or np.ndarray
            Reconstructed signal.
        """
        coef = self.cwt_coeffs[scale]
        peaks = self.find_sequence('YAZB', scale)
        y = np.zeros_like(self.signal)

        for peak in peaks:
            pos = peak * scale
            width = scale
            height = coef[peak]
            y += height * np.exp(-((np.arange(len(y)) - pos) / (0.6006 * width)) ** 2)

        if return_signal:
            from signal import signal as Signal
            return Signal(x=self.x.copy(), y=y, name=f"reconstructed_{self.name}",
                          x_label="x", y_label="reconstructed", y_unit="a.u.")
        return y

    def align_with(self, other, scale=1):
        """
        Align symbolic sequences and compute mutual entropy.

        Returns:
            SimpleNamespace: with fields
                - seq1_aligned (str)
                - seq2_aligned (str)
                - aligned_signal (list of tuples)
                - mutual_entropy (float)
        """
        seq1 = self.codes[scale]['letters']
        seq2 = other.codes[scale]['letters']
        aligner = PairwiseAligner()
        aligner.mode = 'global'
        alignment = aligner.align(seq1, seq2)[0]
        a1 = self.reconstruct_aligned_string(seq1, alignment.aligned[0])
        a2 = self.reconstruct_aligned_string(seq2, alignment.aligned[1])
        pairs = list(zip(a1, a2))
        shared = ''.join(a if a == b else '_' for a, b in pairs)
        H = self.entropy_from_string
        return SimpleNamespace(
            seq1_aligned=a1,
            seq2_aligned=a2,
            aligned_signal=pairs,
            mutual_entropy=H(a1) + H(a2) - 2 * H(shared)
        )

    @staticmethod
    def reconstruct_aligned_string(seq, aligned):
        """Fast reconstruction of aligned signals"""
        result, last = [], 0
        for start, end in aligned:
            result.extend(['-'] * (start - last))
            result.extend(seq[start:end])
            last = end
        return ''.join(result)

    @staticmethod
    def entropy_from_string(s):
        """return the entropy of a string"""
        _, counts = np.unique(list(s), return_counts=True)
        return entropy(counts, base=2)

    @staticmethod
    def print_alignment(seq1, seq2, width=80):
        """print aligned sequences"""
        for i in range(0, len(seq1), width):
            s1 = seq1[i:i+width]
            s2 = seq2[i:i+width]
            match = ''.join('|' if a == b else ' ' for a, b in zip(s1, s2))
            print(s1)
            print(match)
            print(s2)
            print()

    def __repr__(self):
        total = sum(len(c['letters']) for c in self.codes.values()) or 1
        ratio = len(self.signal) / total
        print(f"<DNAsignal(length={len(self.signal)}, scales={self.scales}, "
                f"transforms={len(self.cwt_coeffs)}, compression_ratio={ratio:.2f}, dtype={self.dtype.__name__})>",sep="\n")
        return str(self)

    def __str__(self):
        return f"<DNAsignal(signal_length={len(self.signal)}, filtered={'Yes' if self.filtered_signal is not None else 'No'})>"

    @staticmethod
    def synthetic_signal(x, peaks, baseline=None):
        """Generate flexible synthetic signals. (obsolete)"""
        y = np.zeros_like(x)
        for pos, width, height in peaks:
            y += height * np.exp(-((x - pos) / (0.6006 * width)) ** 2)
        if baseline:
            y += baseline(x)
        return y

    def plot_signals(self, scales=None):
        """Plot signals."""
        if scales is None:
            scales = self.scales
        plt.figure(figsize=(10, 6))
        plt.plot(self.signal, label='Original Signal', linewidth=4, color="k")
        for scale in scales:
            plt.plot(self.cwt_coeffs[scale], label=f'Scale {scale}', linewidth=2, alpha=0.7)
        plt.legend()
        plt.title("Signal and Transformed Scales")
        plt.x_label(f"{self.x_label} [{self.x_unit}]")
        plt.y_label(f"{self.y_label} [{self.y_unit}]")
        plt.show()

    def plot_transforms(self, indices=None, **kwargs):
        """
        Plot the stored CWT-transformed signals as a signal collection.

        Parameters
        ----------
        indices : list[int or str], optional
            Specific scales or names to plot.
        kwargs : passed to `signal_collection.plot`
        """
        if not hasattr(self, "transforms"):
            raise AttributeError("No transformed signals found. Run `compute_cwt()` first.")
        self.transforms.plot(indices=indices, title=f"CWT Transforms: {self.name}", **kwargs)

    def plot_codes(self, scale, ax=None, colormap=None, alpha=0.4):
        """
        Plot the symbolic DNA-like encoding as colored triangle segments.

        Parameters
        ----------
        scale : int
            The scale at which the signal was encoded.
        ax : matplotlib.axes.Axes, optional
            Axis to draw on. If None, a new figure is created.
        colormap : dict, optional
            Custom mapping of letters to colors. Default uses 7 distinct colors.
        alpha : float
            Transparency for the patches. Default is 0.4.
        """

        if ax is None:
            fig, ax = plt.subplots(figsize=(12, 5))
        else:
            fig = plt.gcf()
            ax = plt.gca()

        # Default color mapping
        default_colormap = {
            'A': 'DeepPink','B': 'OrangeRed', 'C': 'Gold',
            'X': 'DeepSkyBlue','Y': 'DodgerBlue', 'Z': 'Purple',
            '_': 'gray'
        }
        cmap = colormap or default_colormap

        # Flags
        cross_flags = self._is_letter_crossing

        # Retrieve segment data
        code_data = self.codes[scale]
        letters = code_data['letters']
        widths = code_data['widths']
        heights = code_data['heights']
        xloc = code_data['xloc']
        coef = self.cwt_coeffs[scale]
        x = self.x if hasattr(self, 'x') and self.x is not None else np.arange(len(coef))

        # Generate path data
        last_y = 0.0
        patchxy = []

        # Process segments with continuity enforcement
        nletters = len(letters)
        for i in range(nletters):
            letter = letters[i]
            repeatedletter_withprevious = letter if i>0 else False
            repeatedletter_withnext = letter if i<(nletters-1) else False
            repeated = repeatedletter_withprevious or repeatedletter_withnext
            y0 = last_y if cross_flags[letter] or repeated else 0.0
            x0, w, h = xloc[i][0], widths[i], heights[i]
            last_y = y0 + h
            x,y = self._get_triangle_from_letter(letter,x0,y0,w,h)
            triangle = [[x[0], y[0]], [x[1], y[1]], [x[2], y[2]], [x[0], y[0]]]
            patchxy.append(triangle)

        # Plot patches
        for letter,verts in zip(letters,patchxy):
            poly = Polygon(verts, closed=True, color=cmap.get(letter, 'black'), alpha=alpha)
            ax.add_patch(poly)

        # Plot raw and transformed signals
        ax.plot(self.x, self.signal, color='black', linewidth=3, label='Original Signal')
        ax.plot(self.x, coef, color='blue', linewidth=2, label=f'CWT Scale {scale}')
        ax.set_xlabel(f"{self.x_label} [{self.x_unit}]")
        ax.set_ylabel(f"{self.y_label} [{self.y_unit}]")
        ax.set_title(f"Symbolic Segments (Scale {scale})")
        ax.legend()
        ax.set_xlim([self.x[0], self.x[-1]])
        ax.axhline(0, color='k', linewidth=0.5, linestyle=':')
        plt.tight_layout()
        plt.show()
        return fig

    def plot_scalogram(self):
        """
        Plot a scalogram with two subplots:
        - Top: colored image of CWT coefficient amplitudes
        - Bottom: line curves of selected scales

        Returns
        -------
        fig : matplotlib.figure.Figure
            The matplotlib figure object.
        """
        X = np.array(self.x)
        scales = list(self.cwt_coeffs.keys())
        coeffs = np.array([self.cwt_coeffs[scale] for scale in scales])

        fig, axs = plt.subplots(2, 1, figsize=(10, 6), gridspec_kw={'height_ratios': [3, 2]}, sharex=True)

        # Top scalogram (image)
        im = axs[0].imshow(coeffs, aspect='auto', extent=[X[0], X[-1], scales[-1], scales[0]],
                           cmap='viridis', interpolation='nearest')
        axs[0].set_ylabel("Scale")
        fig.colorbar(im, ax=axs[0], label="Amplitude")
        axs[0].set_title("Scalogram")

        # Bottom line plots
        for i, scale in enumerate(scales):
            axs[1].plot(X, coeffs[i] + i, label=f"Scale {scale}")
        axs[1].set_xlabel("x")
        axs[1].set_ylabel("Wavelet coefficients")
        axs[1].set_title("Wavelet Coefficients by Scale")
        axs[1].legend(loc="upper right", fontsize='small')

        plt.tight_layout()
        return fig

    @staticmethod
    def _pairwiseEntropyDistance(list_DNAsignals, scale=None,
                  engine=None, engineOpts=None,):
        """
        Calculate excess-entropy pairwise distances.

        Parameters
        ----------
        list_DNAsignals : list of valid DNAsignals (mandatory)
        scale           : int (mandatory)
        engine          : {'difflib', 'bio'} or None (default)
        engineOpts      : dict, optional (alignment parameters for the selected engine)

        Returns
        -------
        D      : nxn np.array of excess entropy distances (H(A) + H(B) - 2 * H(A*B))
        names  : list of signal names
        """

        if not isinstance(list_DNAsignals, list):
            raise TypeError(f"list_DNAsignals must be list not a {type(list_DNAsignals).__name__}")
        for o in list_DNAsignals:
            if not isinstance(o, DNAsignal):
                raise TypeError(f"all listed elements must be a DNAsignal not a {type(o).__name__}")
        if scale is None or not isinstance(scale, int):
            raise TypeError(f"scale must be an int not a {type(scale).__name__}")
        if scale <= 0:
            raise ValueError("scale must be positive")

        n = len(list_DNAsignals)
        D = np.zeros((n, n), dtype=np.float64)

        total_pairs = n * (n - 1) // 2
        pair_index = 0
        start_time = time()

        with tqdm(total=total_pairs, desc="Pairwise distances", unit="pair") as pbar:
            for i in range(n):
                A = list_DNAsignals[i].codesfull[scale]
                for j in range(i):
                    B = list_DNAsignals[j].codesfull[scale]
                    A.align(B, engine=engine, engineOpts=engineOpts)
                    D[i, j] = A.excess_entropy(B)
                    pair_index += 1
                    pbar.update(1)
                    # Optional: show elapsed/ETA in tqdm (already included by default)

        D += D.T
        names = [o.name for o in list_DNAsignals]
        elapsed = time() - start_time
        print(f"Pairwise distances computation completed in {elapsed:.2f} seconds.")
        return DNApairwiseAnalysis(D,names,list_DNAsignals)

    @staticmethod
    def _pairwiseJaccardMotifDistance(list_DNAsignals, scale=None,
                                       pattern='YAZB', minlen=4,
                                       classification='any',  # 'canonical', 'variant', or 'any'
                                       plot=True):
        """
        Compute pairwise Jaccard distances based on motif presence across symbolic DNAstr sequences.

        Parameters
        ----------
        list_DNAsignals : list of DNAsignal
            List of signals with .codesfull at given scale.
        scale : int
            Scale index for accessing .codesfull[scale].
        pattern : str
            Motif pattern to look for (default 'YAZB').
        minlen : int
            Minimum length of valid motifs.
        classification : {'canonical', 'variant', 'any'}
            Filter for motif type.
        plot : bool
            Whether to plot motif positions for each sequence.

        Returns
        -------
        D : np.ndarray
            Symmetric pairwise Jaccard distance matrix.
        names : list of str
            List of DNAsignal names.
        """
        if not isinstance(list_DNAsignals, list):
            raise TypeError("list_DNAsignals must be a list.")
        if scale is None or not isinstance(scale, int) or scale <= 0:
            raise ValueError("scale must be a positive integer.")

        n = len(list_DNAsignals)
        motif_sets = []

        for i, obj in enumerate(list_DNAsignals):
            if not isinstance(obj, DNAsignal):
                raise TypeError(f"Element {i} is not a DNAsignal.")
            dna = obj.codesfull[scale]
            df = dna.extract_motifs(pattern=pattern, minlen=minlen, plot=False)
            if classification == 'canonical':
                df = df[df['classification'] == 'canonical']
            elif classification == 'variant':
                df = df[df['classification'] == 'variant']
            # Represent a sequence as set of motif *start* positions
            motif_set = set(df['start']) if not df.empty else set()
            motif_sets.append(motif_set)

        if plot:
            # Assume all DNAstr have same length, dx, and iloc
            L = len(list_DNAsignals[0].codesfull[scale])
            dx = list_DNAsignals[0].codesfull[scale].dx
            iloc = list_DNAsignals[0].codesfull[scale].iloc
            x_start = iloc * dx if isinstance(iloc, int) else iloc[0] * dx
            x_vals = np.arange(L) * dx + x_start
            prevalence = np.zeros(L, dtype=int)

            for obj in list_DNAsignals:
                dna = obj.codesfull[scale]
                df = dna.extract_motifs(pattern=pattern, minlen=minlen, plot=False)
                if classification == 'canonical':
                    df = df[df['classification'] == 'canonical']
                elif classification == 'variant':
                    df = df[df['classification'] == 'variant']
                for _, row in df.iterrows():
                    start = row['start']
                    width = row['end'] - row['start']  # assuming inclusive-exclusive
                    prevalence[start:start + width] += 1

            plt.figure(figsize=(12, 4))
            plt.plot(x_vals, prevalence/n, marker='o', linestyle='-', alpha=0.8)
            plt.title(f"Motif prevalence at each position for pattern '{pattern}'")
            plt.xlabel("x position")
            plt.ylabel("Number of sequences with motif coverage")
            plt.grid(True)
            plt.tight_layout()
            plt.show()

        D = np.zeros((n, n), dtype=np.float64)
        total_pairs = n * (n - 1) // 2
        pair_index = 0
        start_time = time()

        with tqdm(total=total_pairs, desc="Pairwise Jaccard (motif)", unit="pair") as pbar:
            for i in range(n):
                A = motif_sets[i]
                for j in range(i):
                    B = motif_sets[j]
                    intersection = len(A & B)
                    union = len(A | B)
                    distance = 1.0 if union == 0 else 1 - intersection / union
                    D[i, j] = distance
                    pbar.update(1)

        D += D.T
        names = [o.name for o in list_DNAsignals]
        elapsed = time() - start_time
        print(f"Jaccard motif-based distance computation completed in {elapsed:.2f} seconds.")
        return DNApairwiseAnalysis(D, names, list_DNAsignals)

    @staticmethod
    def _pairwiseJensenShannonDistance(list_DNAsignals, scale=None):
        """
        Calculate pairwise Jensen-Shannon distances between DNAstr codes at a given scale.

        Parameters
        ----------
        list_DNAsignals : list of DNAsignal
            List of valid DNAsignal instances.
        scale : int
            Scale index to select the code from `codesfull`.

        Returns
        -------
        DNApairwiseAnalysis
            Matrix object containing pairwise Jensen-Shannon distances.
        """
        if not isinstance(list_DNAsignals, list):
            raise TypeError("list_DNAsignals must be a list")
        for o in list_DNAsignals:
            if not isinstance(o, DNAsignal):
                raise TypeError(f"All elements must be DNAsignal, not {type(o).__name__}")
        if scale is None or not isinstance(scale, int):
            raise TypeError("scale must be an integer")
        if scale < 0:
            raise ValueError("scale must be non-negative")

        n = len(list_DNAsignals)
        D = np.zeros((n, n), dtype=np.float64)
        total_pairs = n * (n - 1) // 2
        start_time = time()

        with tqdm(total=total_pairs, desc="Jensen-Shannon", unit="pair") as pbar:
            for i in range(n):
                A = list_DNAsignals[i].codesfull[scale]
                for j in range(i):
                    B = list_DNAsignals[j].codesfull[scale]
                    D[i, j] = A.jensen_shannon(B)
                    pbar.update(1)

        D += D.T
        names = [o.name for o in list_DNAsignals]
        elapsed = time() - start_time
        print(f"Jensen-Shannon distance matrix completed in {elapsed:.2f} seconds.")
        return DNApairwiseAnalysis(D, names, list_DNAsignals)

    @staticmethod
    def _pairwiseLevenshteinDistance(list_DNAsignals, scale=None,
                                      use_alignment=False,
                                      engine=None,
                                      engineOpts=None,
                                      forced=False):
        """
        Compute pairwise Levenshtein distances between codes at a given scale.

        Parameters
        ----------
        list_DNAsignals : list of DNAsignal
            List of DNAsignal objects to compare.
        scale : int
            Scale index in codesfull.
        use_alignment : bool, optional
            If True, align codes before computing distance. Default is Full.
        engine : str, optional
            Alignment engine ('difflib' or 'bio') if use_alignment is True.
        engineOpts : dict, optional
            Parameters for the alignment engine.
        forced : bool, optional
            If True, allow forced alignment even if dx mismatches.

        Returns
        -------
        DNApairwiseAnalysis
            Object holding the nxn Levenshtein distance matrix.
        """
        if not isinstance(list_DNAsignals, list):
            raise TypeError("list_DNAsignals must be a list")
        for o in list_DNAsignals:
            if not isinstance(o, DNAsignal):
                raise TypeError(f"All elements must be DNAsignal, not {type(o).__name__}")
        if scale is None or not isinstance(scale, int):
            raise TypeError("scale must be an integer")
        if scale < 0:
            raise ValueError("scale must be non-negative")

        n = len(list_DNAsignals)
        D = np.zeros((n, n), dtype=np.float64)
        total_pairs = n * (n - 1) // 2
        start_time = time()

        with tqdm(total=total_pairs, desc="Levenshtein", unit="pair") as pbar:
            for i in range(n):
                A = list_DNAsignals[i].codesfull[scale]
                for j in range(i):
                    B = list_DNAsignals[j].codesfull[scale]
                    D[i, j] = A.levenshtein(B,
                                            use_alignment=use_alignment,
                                            engine=engine,
                                            engineOpts=engineOpts,
                                            forced=forced)
                    pbar.update(1)

        D += D.T
        names = [o.name for o in list_DNAsignals]
        elapsed = time() - start_time
        print(f"Levenshtein distance matrix completed in {elapsed:.2f} seconds.")
        return DNApairwiseAnalysis(D, names, list_DNAsignals)


    @property
    def letters(self):
        """Return used letters"""
        if hasattr(self,"sincodesfull") and isinstance(self.sincodesfull,DNAFullCodes) and len(self.scales)>0:
            return list(self.sincodesfull[self.scales[0]].keys())
        elif hasattr(self,"codesfull") and isinstance(self.codesfull,DNAFullCodes) and len(self.scales)>0:
            return list(set(str(self.codesfull[self.scales[0]])))
        elif hasattr(self,"codes") and isinstance(self.codes,DNACodes) and len(self.scales)>0:
            return list(set(self.codes[self.scales[0]]["letters"]))
        else:
            return None

# ------------------------
# DNAstr class
# ------------------------
class DNAstr(str):
    """
    A symbolic DNA-like sequence class supporting alignment, entropy analysis,
    edit-distance metrics, and signal reconstruction from symbolic codes.

    Extended from `str`, it is designed for symbolic transformations of signals
    (e.g., wavelet-encoded GC-MS peaks or time series).

    Main Features
    -------------
    - Supports symbolic operations for pattern recognition, entropy, alignment.
    - Encodes x-resolution (`dx`), original index (`iloc`), and physical x-range (`xloc`).
    - Aligns sequences with visual inspection and rich diffs.
    - Converts symbolic strings into synthetic numerical signals.

    Operators
    ---------
    + : concatenate two DNAstr objects
    - : symbolic difference after alignment (mismatches only)
    == : equality comparison (exact content and dx)

    Key Methods
    -----------
    - __init__ / __new__      : Constructor with metadata (`dx`, `iloc`, `xloc`)
    - align(other)            : Align this DNAstr to another, update mask and aligned views
    - wrapped_alignment()     : Pretty terminal view of the alignment with colors and symbols
    - html_alignment()        : Rich HTML display of the alignment (Jupyter)
    - plot_alignment()        : Visualize waveform alignment with symbolic signals
    - plot_mask()             : Color block plot showing matches/mismatches/gaps
    - find(pattern, regex=False) : Search for symbolic patterns with fuzziness or regex
    - to_signal()             : Convert symbolic code into synthetic signal (NumPy)
    - vectorized()            : Convert string to integer codes
    - summary()               : Print entropy and character frequencies
    - mutation_counts         : Property: {'matches', 'mismatches', 'indels'}
    - entropy                 : Property: Shannon entropy
    - mutual_entropy(other)   : Mutual entropy of two sequences
    - excess_entropy(other)   : Excess entropy H1 + H2 - 2 * H12
    - jensen_shannon(other)   : Jensen-Shannon divergence
    - jaccard(other)          : Jaccard similarity
    - alignment_stats         : Property: Match, substitution, gap counts
    - score(normalized=True) : Alignment score (fraction of matches)
    - has(other: str)         : Check if a pattern or substring exists

    Attributes
    ----------
    dx : float
        Average resolution along the x-axis.
    iloc : int or tuple of int
        Positional index or index range in the source DNA string.
    xloc : float or tuple of float
        Corresponding x-value(s) for the symbolic sequence.
    aligned_with : str or None
        Aligned form of self with insertions (spaces) where needed.
    other_copy : str or None
        Aligned form of the reference sequence.
    ref_hash : str or None
        SHA256 hash of the aligned reference sequence.
    mask : str or None
        Alignment mask: '=' for matches, '*' for substitutions, ' ' for gaps.
    engine : str
        Alignment engine: 'difflib' or 'bio'.
    engineOpts : dict
        Options passed to the alignment engine.

    Examples
    --------
    >>> s1 = DNAstr("YYAAZZBB", dx=0.5)
    >>> s2 = DNAstr("YAABZBB", dx=0.5)
    >>> s1.align(s2)
    >>> print(s1.wrapped_alignment(40))
    >>> s1.plot_alignment()
    >>> segments = s1.find("YAZB")
    >>> segments[0].to_signal().plot()

    """

    def __new__(cls, content, dx=1.0, iloc=0, xloc=None, x_label="index", x_unit="-",
                engine="difflib", engineOpts=None):
        """
        Construct a new DNAstr object.

        Parameters
        ----------
        content : str
            The symbolic DNA-like string content.
        dx : float, optional
            Nominal resolution (default: 1.0).
        iloc : int, optional
            Integer index or start index of the sequence (default: 0).
        xloc : float, optional
            X-coordinate of the sequence origin.
        engine : {'difflib', 'bio'}, optional
            Default alignment engine to use.
        engineOpts : dict, optional
            Dictionary of alignment parameters for the selected engine.

        Returns
        -------
        DNAstr
            Initialized DNAstr instance.
        """
        obj = str.__new__(cls, content)
        obj.dx = dx
        obj.iloc = iloc
        obj.xloc = xloc
        obj.x_label = x_label
        obj.x_unit = x_unit
        obj.aligned_with = None
        obj.other_copy = None
        obj.ref_aligned = None
        obj.ref_hash = None
        obj.mask = None
        if engine not in ("difflib", "bio"):
            raise ValueError('engine must be "difflib" or "bio"')
        obj.engine = engine
        obj.engineOpts = engineOpts or {}
        obj._hash = hashlib.sha256(obj.encode()).hexdigest()
        return obj

    def __hash__(self):
        """
        Return hash combining the string content and dx.

        Returns
        -------
        int
            Hash of the DNAstr object.
        """
        return hash((str(self), self.dx))

    def summary(self):
        """
        Summarize the DNAstr with key stats: length, unique letters, entropy, etc.

        Returns
        -------
        dict
            Dictionary containing length, letter frequency, Shannon entropy, and dx.
        """
        length = len(self)
        freqs = Counter(self)
        prob = np.array(list(freqs.values())) / length
        entropy = -np.sum(prob * np.log2(prob))
        return {
            'length': length,
            'letters': dict(freqs),
            'entropy (Shannon)': entropy,
            'dx': self.dx
        }

    def vectorized(self, codebook={"A":1,"B":2,"C":3,"X":4,"Y":5,"Z":6,"_":0}):
        """
        Map the DNAstr content to an integer array using a codebook.

        Parameters
        ----------
        codebook : dict, optional
            Dictionary mapping characters to integer values.
            default = {"A":1,"B":2,"C":3,"X":4,"Y":5,"Z":6,"_":0}
            None will generate a codebook based on current symbols only

        Returns
        -------
        np.ndarray
            Vectorized integer representation of the string.
        """
        from collections import OrderedDict
        if codebook is None:
            unique_chars = list(OrderedDict.fromkeys(self))
            codebook = {c: i for i, c in enumerate(unique_chars)}
        return np.array([codebook.get(c, -1) for c in self], dtype=int)

    def __eq__(self, other):
        """
        Check equality based on symbolic content and dx resolution.

        Parameters
        ----------
        other : DNAstr
            Another DNAstr instance.

        Returns
        -------
        bool
            True if both content and dx are equal.
        """
        return isinstance(other, DNAstr) and str.__eq__(self, other) and self.dx == other.dx


    def __add__(self, other):
        """
        Concatenate two DNAstr instances with identical dx values.

        Parameters
        ----------
        other : DNAstr
            Another DNAstr sequence.

        Returns
        -------
        DNAstr
            Concatenated DNAstr sequence.

        Raises
        ------
        TypeError
            If the argument is not a DNAstr.
        ValueError
            If `dx` values differ.
        """
        if not isinstance(other, DNAstr):
            raise TypeError("Can only concatenate DNAstr (not '{}')".format(type(other).__name__))
        if self.dx != other.dx:
            raise ValueError("dx mismatch. Use forced=True to override.")
        return DNAstr(str.__add__(self, other), dx=self.dx, iloc=self.iloc, xloc=self.xloc,
                                   x_label=self.x_label, x_unit=self.x_unit)

    def __sub__(self, other):
        """
        Subtract two DNAstr sequences by aligning and removing matched regions.

        Parameters
        ----------
        other : DNAstr
            DNAstr instance to align and subtract.

        Returns
        -------
        DNAstr
            A new DNAstr containing mismatched symbols only.
        """
        if not isinstance(other, DNAstr):
            raise TypeError("Can only subtract DNAstr instances")
        self.align(other)
        mismatch = ''.join([a for a, b in zip(self.aligned_with, self.other_copy) if a != b and b != ' '])
        return DNAstr(mismatch, dx=self.dx, iloc=self.iloc, xloc=self.xloc,
                                   x_label=self.x_label, x_unit=self.x_unit)

    @property
    def mutation_counts(self):
        """Counts of insertions, deletions/substitutions, and matches."""
        m = self.mask
        return {
            'matches': m.count('='),
            'mismatches': m.count('*'),
            'indels': m.count(' ')
        }
    @property
    def entropy(self):
        """Compute the Shannon entropy of the DNAstr sequence"""
        count = Counter(self)
        total = sum(count.values())
        return -sum((v / total) * np.log2(v / total) for v in count.values())

    def mutual_entropy(self, other=None):
        """Compute the Shannon mutual entropy of two DNAstr sequences from their aligned segments"""
        if other is not None and not isinstance(other,DNAstr):
            raise TypeError(f"other must be a DNAstr not a {type(self).__name__}")
        if other is None and (not hasattr(self,"aligned_with") or self.aligned_with is None):
            raise ValueError("align the code with .align(other) or provide other")
        if isinstance(other,DNAstr):
            self.align(other)
        aligned = self.aligned_code
        count = Counter(aligned)
        total = sum(count.values())
        return -sum((v / total) * np.log2(v / total) for v in count.values())

    def excess_entropy(self, other):
        """Compute the excess Shannon entropy of two DNAstr sequences H(A)+H(B)-2*H(AB)"""
        return self.entropy + other.entropy - 2 * self.mutual_entropy(other)

    def jensen_shannon(self, other, base=2):
        """
        Compute the Jensen-Shannon distance between self and another DNAstr.

        Parameters
        ----------
        other : DNAstr
            Another DNAstr instance.
        base : float, optional
            Base for the logarithm (default: 2)

        Returns
        -------
        float
            Jensen-Shannon distance.
        """
        v1 = Counter(self)
        v2 = Counter(other)
        all_keys = sorted(set(v1) | set(v2))
        p = np.array([v1.get(k, 0) for k in all_keys], dtype=float)
        q = np.array([v2.get(k, 0) for k in all_keys], dtype=float)
        p /= p.sum()
        q /= q.sum()
        return jensenshannon(p, q, base=base)


    def levenshtein(self, other, use_alignment=True, engine=None, engineOpts=None, forced=False):
        """
        Compute the Levenshtein distance between this DNAstr and another one.

        Parameters
        ----------
        other : DNAstr
            Another DNAstr object to compare against.
        use_alignment : bool, default=True
            If True, uses the aligned sequences (computed if necessary).
            If False, compares the raw sequences directly.
        engine : {'difflib', 'bio'}, optional
            Alignment engine to use if alignment is needed.
        engineOpts : dict, optional
            Parameters for the selected alignment engine.
        forced : bool, default=False
            Force alignment even if dx values differ.

        Returns
        -------
        dist : int
            Levenshtein distance between the two sequences (aligned or raw).

        Examples
        --------
        A = DNAstr("YAZBZAY")
        B = DNAstr("YAZBZZY")
        A.levenshtein_distance(B, use_alignment=False)  # raw
        A.levenshtein_distance(B, use_alignment=True, engine="bio")  # aligned
        """
        if not isinstance(other, DNAstr):
            raise TypeError("Argument must be a DNAstr instance")
        if use_alignment:
            self.align(other, engine=engine, engineOpts=engineOpts, forced=forced)
            s1, s2 = self.aligned_with, self.other_copy
        else:
            s1, s2 = str(self), str(other)
        return Levenshtein.distance(s1, s2)


    def jaccard(self, other):
        """
        Compute the Jaccard distance between two DNAstr sequences.

        Parameters
        ----------
        other : DNAstr
            The other DNAstr sequence to compare with.

        Returns
        -------
        float
            Jaccard distance: 1 - (intersection / union) of unique letters.
        """
        set_self = set(self)
        set_other = set(other)
        intersection = set_self & set_other
        union = set_self | set_other
        return 1 - len(intersection) / len(union) if union else 0.0

    def align(self, other, engine=None, engineOpts=None, forced=False):
        """
        Align this DNAstr sequence to another, allowing insertions/deletions to maximize matches.

        Parameters
        ----------
        other : DNAstr
            Another DNAstr object to align with.
        engine : {'difflib', 'bio'} or None
            Alignment engine to use:
                - 'difflib': uses difflib.SequenceMatcher (fast, approximate).
                - 'bio'   : uses Bio.Align.PairwiseAligner (biologically inspired global alignment).
            If None, defaults to self.engine.
        engineOpts : dict, optional
            Dictionary of alignment parameters for the selected engine.
        forced : bool
            If True, allow alignment even if `dx` values differ. If False (default), a mismatch in
            `dx` will raise an error to prevent incorrect alignment of signals with different sampling.

        Returns
        -------
        aligned_self : str
            Aligned version of this sequence (with gaps inserted where needed).
        aligned_other : str
            Aligned version of the other sequence.

        Notes
        -----
        The alignment is symmetric and permanent: both sequences are aligned with
        gaps introduced (spaces) to preserve positional correspondence. A hash of
        the aligned `other` sequence is stored to detect redundant alignments.

        A match mask (`self.mask`) is generated with:
            '=' for exact matches,
            '*' for mismatches (substitutions),
            ' ' for insertions/deletions (gaps).

        The method updates:
            - self.aligned_with
            - self.other_copy
            - self.mask
            - self.ref_hash

        Example:
        --------
        S1 = DNAstr("AABBCC")
        S2 = DNAstr("AACBCC")
        S1.align(S2,"difflib")
        print(S1.mask)
        print(S1.wrapped_alignment())
        ==*===
        AACBCC
        || |||
        AABBCC

        S1 = DNAstr("AABBCC")
        S2 = DNAstr("AACBCC")
        S1.align(S2,"bio")
        print(S1.mask)
        print(S1.wrapped_alignment())
        ==  ==
        AAB·CC
        ||  ||
        AA·BCC

        S1 = DNAstr("AABBCCXYZZZ")
        S2 = DNAstr("AACBCCZZXXX")
        S1.align(S2,"bio")
        print(S1.mask)
        print(S1.wrapped_alignment())
        == *   ==
        AABCC··ZZ
        ||     ||
        AA·B·CCZZ

        """
        if not isinstance(other, DNAstr):
            raise TypeError("Alignment requires another DNAstr instance")
        if not forced and self.dx != other.dx:
            raise ValueError("dx mismatch. Use forced=True to override.")
        if hasattr(other, 'ref_hash') and self.ref_hash is not None and self.ref_hash == other.ref_hash:
            return self.aligned_with, self.other_copy

        engine = engine or self.engine
        engineOpts = engineOpts or self.engineOpts.get(engine, {})

        if engine == 'difflib':
            sm = SequenceMatcher(None, other, self)
            aligned_self, aligned_other = [], []
            for tag, i1, i2, j1, j2 in sm.get_opcodes():
                if tag == 'equal':
                    aligned_self.extend(self[j1:j2])
                    aligned_other.extend(other[i1:i2])
                elif tag == 'replace':
                    aligned_self.extend(self[j1:j2])
                    aligned_other.extend(other[i1:i2])
                elif tag == 'insert':
                    aligned_self.extend(self[j1:j2])
                    aligned_other.extend(' ' * (j2 - j1))
                elif tag == 'delete':
                    aligned_self.extend(' ' * (i2 - i1))
                    aligned_other.extend(other[i1:i2])

        elif engine == 'bio':
            aligner = PairwiseAligner()
            for k, v in engineOpts.items():
                setattr(aligner, k, v)

            alignment = aligner.align(other, self)[0]  # Best alignment
            aligned_self = []
            aligned_other = []

            # These are lists of (start, end) index tuples for each sequence
            self_blocks = alignment.aligned[1]
            other_blocks = alignment.aligned[0]

            self_pos = 0
            other_pos = 0

            for (o_start, o_end), (s_start, s_end) in zip(other_blocks, self_blocks):
                # Fill gaps in other
                if o_start > other_pos:
                    gap_len = o_start - other_pos
                    aligned_other.extend(other[other_pos:o_start])
                    aligned_self.extend([' '] * gap_len)
                    other_pos = o_start
                # Fill gaps in self
                if s_start > self_pos:
                    gap_len = s_start - self_pos
                    aligned_self.extend(self[self_pos:s_start])
                    aligned_other.extend([' '] * gap_len)
                    self_pos = s_start

                # Aligned regions
                aligned_self.extend(self[s_start:s_end])
                aligned_other.extend(other[o_start:o_end])
                self_pos = s_end
                other_pos = o_end

            # Tail padding
            aligned_self.extend(self[self_pos:])
            aligned_other.extend([' '] * (len(self) - self_pos))
            aligned_other.extend(other[other_pos:])
            aligned_self.extend([' '] * (len(other) - other_pos))

        else:
            raise ValueError("Unknown alignment engine: choose 'difflib' or 'bio'")

        self.aligned_with = ''.join(aligned_self)
        self.other_copy = ''.join(aligned_other)
        if len(self.aligned_with) != len(self.other_copy):
            raise RuntimeError("Mismatch in alignment lengths: check alignment logic.")
        self.mask = ''.join('=' if a == b else '*' if b != ' ' and a != ' ' else ' '
                            for a, b in zip(self.aligned_with, self.other_copy))
        self.ref_hash = hashlib.sha256(self.other_copy.encode()).hexdigest()
        self.engine = engine
        self.engineOpts[engine] = engineOpts
        return self.aligned_with, self.other_copy

    @property
    def alignment_stats(self):
        """Retrun DNAstr alignment statistics"""
        if self.mask is None:
            raise ValueError("No alignment performed yet.")
        return {
            "matches": self.mask.count('='),
            "substitutions": self.mask.count('*'),
            "gaps": self.mask.count(' ')
        }

    @property
    def aligned_code(self):
        """return aligned code"""
        if not hasattr(self,"aligned_with") or self.aligned_with is None:
            raise ValueError("the code is not aligned")
        return re.sub(r'[^A-CX-Z]', '', self.aligned_with)

    def score(self, normalized=True):
        """
        Return an alignment score, optionally normalized.

        Parameters
        ----------
        normalized : bool
            If True (default), return score as a fraction of total aligned positions.

        Returns
        -------
        float
            Alignment score.
        """
        stats = self.alignment_stats
        score = stats["matches"]
        return score / len(self.mask) if normalized else score

    @staticmethod
    def _supports_color():
        """Returns True if ther terminal supports colors"""
        return hasattr(sys.stdout, "isatty") and sys.stdout.isatty() and os.getenv("TERM") not in (None, "dumb")

    def wrapped_alignment(self, width=80, colors=True):
        """
        Return a line-wrapped alignment view (multi-line), optionally color-coded
        for terminal/IPython usage (Spyder, Jupyter).

        Parameters
        ----------
        width : int
            Number of characters per line in wrapped display.
        colors : bool
            If True, use ANSI codes to highlight differences. May be overridden
            if terminal does not support ANSI (e.g., Spyder).

        Returns
        -------
        str
            Wrapped, optionally colorized alignment.
        """
        if self.aligned_with is None or self.other_copy is None:
            raise ValueError("Alignment has not been computed yet.")

        match_mask = self.mask
        s1 = self.other_copy
        s2 = self.aligned_with

        if colors and not DNAstr._supports_color():
            colors = False

        def colorize(c, match):
            if not colors:
                return c
            if c == ' ':
                return '\x1b[90m·\x1b[0m'
            elif match == '|':
                return f'\x1b[92m{c}\x1b[0m'
            else:
                return f'\x1b[91m{c}\x1b[0m'

        lines = []
        for i in range(0, len(s1), width):
            s1_block = s1[i:i+width]
            s2_block = s2[i:i+width]
            msk_block = match_mask[i:i+width]
            s1c = ''.join(colorize(c, m) for c, m in zip(s1_block, msk_block))
            s2c = ''.join(colorize(c, m) for c, m in zip(s2_block, msk_block))
            match_line = ''.join('|' if m == '=' else ' ' for m in msk_block)
            lines.extend([s1c, match_line, s2c, ''])
        return '\n'.join(lines)

    def html_alignment(self):
        """
        Render the alignment using HTML with color coding:
        - green: match
        - blue: gap
        - red: substitution

        Returns
        -------
        None
            Displays HTML directly in Jupyter/Notebook environments.
        """
        if not self.aligned_with or not self.other_copy:
            raise ValueError("Alignment not available. Call .align() first.")
        html = "<pre style='font-family: monospace;'>"
        for a, b in zip(self.other_copy, self.aligned_with):
            if a == b:
                html += f"<span style='color:green'>{b}</span>"
            elif a == ' ' or b == ' ':
                html += f"<span style='color:blue'>{b}</span>"
            else:
                html += f"<span style='color:red'>{b}</span>"
        html += "</pre>"
        display(HTML(html))

    def __repr__(self):
        """
        Return a short technical representation of the DNAstr instance.

        Returns
        -------
        str
            Description of alignment status and length.
        """
        base = f"<DNAstr: {len(self)} symbols"
        if self.aligned_with:
            base += f" - aligned against <HASH {self.ref_hash[:8]}>>"
        else:
            base += " - not aligned>"
        return base

    def __str__(self):
        """
        String representation.

        Returns
        -------
        str
            Original string content.
        """
        #return repr(self) if self.aligned_with else str.__str__(self)
        return super().__str__()

    def find(self, pattern, regex=False):
        """
        Finds all fuzzy (or regex-based) occurrences of a DNA-like sequence pattern.

        Parameters
        ----------
        pattern : str
            The symbolic sequence to search for (e.g., "YAZB").
        regex : bool, optional
            If False (default), interprets pattern as symbolic and inserts '.' between characters.
            If True, uses the raw pattern as a regular expression.

        Returns
        -------
        list of DNAstr
            A list of DNAstr slices with attributes:
                - iloc: (start_idx, end_idx)
                - xloc: (x_start, x_end)
                - width: segment width
        """
        if not regex:
            # Turn 'YAZB' into 'Y+A+Z+B+?'
            pattern = ''.join(f"{c}+" for c in pattern)  # Greedy
        matches = []
        for m in re.finditer(pattern, str(self)):
            start, end = m.span()
            substr = self[start:end]
            dna = DNAstr(substr,
                         dx=self.dx,
                         iloc=(start,end),
                         xloc=(self.xloc[0]+start*self.dx,self.xloc[0]+end*self.dx),
                         x_label=self.x_label, x_unit=self.x_unit)
            dna.iloc = (start, end)
            if hasattr(self, "xloc") and self.xloc is not None:
                x_start = self.xloc[0] + self.dx * start
                x_end = self.xloc[0] + self.dx * end
                dna.xloc = (x_start, x_end)
            else:
                dna.xloc = (start * self.dx, end * self.dx)
            dna.width = end - start
            matches.append(dna)
        return matches


    def to_signal(self):
        """
        Converts the symbolic DNA sequence into a synthetic NumPy array mimicking the original wavelet-transformed signal.

        Rules per letter:
            - 'A': Crosses zero upward → linear from -1 to +1, zero in the middle
            - 'Z': Crosses zero downward → linear from +1 to -1, zero in the middle
            - 'B': Increasing negative → from -1 to 0
            - 'Y': Decreasing negative → from 0 to -1
            - 'C': Increasing positive → from 0 to +1
            - 'X': Decreasing positive → from +1 to 0
            - '_': Flat at 0

        Returns
        -------
        numpy.ndarray
            Synthetic signal array matching the symbolic encoding.
        """
        s = []
        i = 0
        while i < len(self):
            letter = self[i]
            j = i
            while j < len(self) and self[j] == letter:
                j += 1
            width = j - i
            if width < 2:
                i = j
                continue  # skip invalid segments
            if letter == 'A':
                seg = np.linspace(-1.0, 1.0, width)
            elif letter == 'Z':
                seg = np.linspace(1.0, -1.0, width)
            elif letter == 'B':
                seg = np.linspace(-1.0, 0, width)
            elif letter == 'Y':
                seg = np.linspace(0, -1.0, width)
            elif letter == 'C':
                seg = np.linspace(0, 1.0, width)
            elif letter == 'X':
                seg = np.linspace(1, 0, width)
            else:  # '_'
                seg = np.zeros(width)

            s.append(seg)
            i = j
        y = np.concatenate(s)
        if self.xloc is None:
            x = None
        else:
            n = len(y)
            x0 = self.xloc[0] if isinstance(self.xloc,(tuple,list)) else self.xloc
            x = np.linspace(x0,x0+self.dx*n,n,endpoint=True,dtype=type(x0))
        return signal(x=x,y=y,name=self._hash)

    def plot_mask(self):
        """
        Plot a color-coded mask of the alignment between sequences.

        Returns
        -------
        matplotlib.figure.Figure
            Matplotlib figure of the alignment mask.
        """
        if not self.other_copy:
            raise ValueError("Alignment required for plotting.")
        fig, ax = plt.subplots(figsize=(12, 2))
        colors = {'=': 'green', '*': 'red', ' ': 'gray'}
        for i, (a, b, m) in enumerate(zip(self.aligned_with, self.other_copy, self.mask)):
            ax.add_patch(Rectangle((i, 0), 1, 1, color=colors[m]))
        ax.set_xlim(0, len(self.aligned_with))
        ax.set_yticks([])
        ax.set_title("DNAstr Alignment Mask")
        ax.set_xlabel("Position")
        return fig

    def plot_alignment(self, dx=1.0, dy=1.0, width=20, normalize=True):
        """
        Plot a block alignment view of two DNAstr sequences with color-coded segments.

        Parameters
        ----------
        dx : float
            Horizontal step between segments (defaults to 1.0).
        dy : float
            Vertical height increment for symbolic waveform visualization.
        width : int
            Number of characters per row (line wrapping).

        Returns
        -------
        matplotlib.figure.Figure
        matplotlib.axes.Axes
        """
        if not self.other_copy:
            raise ValueError("Alignment required for plotting.")
        aligned_self, aligned_other, mask = self.aligned_with, self.other_copy, self.mask
        n = len(aligned_self)
        fig, ax = plt.subplots(figsize=(12, 3))

        def letter_to_height(letter, base=0):
            """Simple deterministic up/down movement from symbolic codes."""
            return base + {
                'A': +1, 'B': +1, 'C': +1,
                'X': -1, 'Y': -1, 'Z': -1,
                '_': 0, ' ': 0
            }.get(letter, 0) * dy

        x, y1, y2 = 0, -0.5, -4
        xs, ys1, ys2 = [0], [0], [y2]

        for i in range(n):
            ax.add_patch(Rectangle((x, y1), dx, dy, color='lightgreen' if mask[i] == '=' else 'salmon', alpha=0.6))
            ys1.append(letter_to_height(aligned_self[i], ys1[-1]))
            ys2.append(letter_to_height(aligned_other[i], ys2[-1]))
            xs.append(x + dx)
            x += dx

        # Convert to numpy arrays and normalize
        if normalize:
            xs = np.array(xs)
            ys1 = np.array(ys1)
            ys2 = np.array(ys2)
            ymax = max(abs(ys1).max(), abs(ys2).max())
            scale = ymax if ymax>1e-12 else 1.0
            ys1 = ys1 / scale
            ys2 = ys2 / scale

        ax.plot(xs, ys1, label="Self", color="DarkMagenta", linewidth=4)
        ax.plot(xs, ys2, label="Reference", color="DodgerBlue", linestyle='-', linewidth=4)

        ax.set_title("Waveform Alignment Visualization")
        ax.set_xlabel("Position")
        ax.set_ylabel("Symbolic Signal")
        ax.legend()
        ax.grid(True)
        plt.tight_layout()
        return fig, ax

    def extract_motifs(self, pattern='YAZB', minlen=4, plot=True):
        """
        Extract and analyze YAZB motifs (canonical and distorted) from the symbolic sequence.

        Parameters
        ----------
        pattern : str
            Canonical motif pattern (default is 'YAZB').
        minlen : int
            Minimum motif length to be considered valid.
        plot : bool
            If True, generate a motif density plot using xloc or sequence index.

        Returns
        -------
        pd.DataFrame
            Table of detected motifs with start/end positions, length, and classification.
        """
        sequence = str(self)
        canonical = pattern
        motif_re = re.compile(r'Y+A+Z+B+')

        matches = []
        for m in motif_re.finditer(sequence):
            start, end = m.span()
            substr = m.group()
            motif_len = end - start
            canonical_match = substr == canonical
            classification = 'canonical' if canonical_match else 'variant'
            if motif_len >= minlen:
                matches.append({
                    'start': start,
                    'end': end,
                    'length': motif_len,
                    'sequence': substr,
                    'classification': classification
                })

        df = pd.DataFrame(matches)

        if plot and not df.empty:
            # Basic position handling
            if hasattr(self, 'xloc') and isinstance(self.xloc, (tuple, list)) and len(self.xloc) == 2:
                x0, x1 = self.xloc
                xspan = np.linspace(x0, x1, len(sequence)+1)
                df['x_start'] = df['start'].apply(lambda i: xspan[i])
                df['x_end'] = df['end'].apply(lambda i: xspan[i])
                xvals = xspan
            else:
                df['x_start'] = df['start']
                df['x_end'] = df['end']
                xvals = np.arange(len(sequence))

            fig, ax = plt.subplots(figsize=(12, 2.5))
            ax.plot(xvals, [1]*len(xvals), alpha=0.1)  # Background for alignment
            for _, row in df.iterrows():
                color = 'green' if row['classification'] == 'canonical' else 'orange'
                ax.axvspan(row['x_start'], row['x_end'], color=color, alpha=0.4)
            ax.set_title(f"Motif Regions in `{getattr(self, '_hash', 'DNAstr')}`")
            ax.set_yticks([])
            ax.set_xlabel("Position (xloc or index)")
            ax.set_xlim([xvals[0], xvals[-1]])
            ax.grid(True, axis='x', linestyle='--', alpha=0.3)
            plt.tight_layout()
            plt.show()
            return df,fig

        return df

# --------------------------
# DNAPairwiseAnalysis class
# --------------------------

class DNApairwiseAnalysis:
    """
    Class to handle pairwise distance analysis, PCoA, clustering, and visualization
    for DNA-coded signals.

    Attributes
    ----------
    D : np.ndarray
        Pairwise excess entropy distance matrix.
    names : list
        Names of the DNA signals.
    DNAsignals : list
        original DNAsignal objects
    coords : np.ndarray
        Coordinates in reduced space (PCoA).
    dimensions : list
        Selected dimensions for reduced analysis.
    linkage_matrix : np.ndarray
        Linkage matrix used for hierarchical clustering.
    """

    def __init__(self, D, names, DNAsignals, name=None):
        self.name = name if not name is None else "unamed"
        self.D = np.array(D)
        self.names = list(names)
        self.n = len(self.names)
        self.DNAsignals = DNAsignals
        self.coords = None
        self.dimensions = list(range(self.n))  # All by default
        self.linkage_matrix = None
        self.pcoa()

    def pcoa(self, n_components=None):
        """Perform Principal Coordinate Analysis (PCoA)."""
        if n_components is None:
            n_components = min(self.n, 1000)
        mds = MDS(n_components=n_components, dissimilarity='precomputed', random_state=42)
        self.coords = mds.fit_transform(self.D)
        self.dimensions = list(range(n_components))

    def select_dimensions(self, dims):
        """Update active dimensions."""
        if isinstance(dims,int):
            self.dimensions = list(range(dims))
        elif isinstance(dims,(list,tuple)):
            self.dimensions = dims
        else:
            raise TypeError(f"dims must be a int, list or tuple not a {type(dims).__name__}")

    def reduced_distances(self):
        """Recompute distances on selected subspace."""
        coords_sel = self.coords[:, self.dimensions]
        return pairwise_distances(coords_sel)

    def compute_linkage(self, method='ward'):
        """Compute hierarchical clustering."""
        D = self.reduced_distances()
        D = (D + D.T) / 2  # force symmetry
        reduced_D = squareform(D)
        self.linkage_matrix = linkage(reduced_D, method=method)

    def heatmap(self, figsize=(10, 8)):
        """Plot heatmap of pairwise distances."""
        fig=plt.figure(figsize=figsize)
        sns.heatmap(self.D, xticklabels=self.names, yticklabels=self.names, cmap='viridis')
        title = f'Pairwise Excess Entropy Distances — {getattr(self, "name", "Unnamed Analysis")}'
        plt.title(title)
        plt.tight_layout()
        plt.show()
        return fig

    def get_cluster_labels(self, n_clusters=2, method='ward'):
        """
        Returns cluster labels from hierarchical clustering. If not computed yet, computes linkage.

        Parameters
        ----------
        n_clusters : int
            Number of clusters to assign.
        method : str
            Linkage method to use if recomputing linkage.

        Returns
        -------
        labels : np.ndarray of int
            Cluster IDs for each sample.
        """
        if self.linkage_matrix is None:
            self.compute_linkage(method=method)
        self.cluster_labels = fcluster(self.linkage_matrix, t=n_clusters, criterion='maxclust')
        return self.cluster_labels

    def scatter(self, dims=(0, 1), annotate=True, figsize=(8, 6), n_clusters=None):
        """
        2D scatter plot in selected dimensions with optional cluster-based coloring.

        Parameters
        ----------
        dims : tuple
            Dimensions to plot (default: (0, 1)).
        annotate : bool
            If True, annotate points with their index.
        figsize : tuple
            Size of the plot.
        n_clusters : int or None
            If provided, use clustering to color points.

        Returns
        -------
        fig : matplotlib.figure.Figure
        """
        fig = plt.figure(figsize=figsize)
        x, y = self.coords[:, dims[0]], self.coords[:, dims[1]]

        if n_clusters is not None:
            labels = self.get_cluster_labels(n_clusters=n_clusters)
            scatter = plt.scatter(x, y, c=labels, cmap='tab10', s=50)
        else:
            scatter = plt.scatter(x, y, s=50)

        if annotate:
            for i, name in enumerate(self.names):
                plt.text(x[i], y[i], str(i), fontsize=9)

        plt.xlabel(f"Dimension {dims[0] + 1}")
        plt.ylabel(f"Dimension {dims[1] + 1}")
        title = f'PCoA Scatter Plot — {getattr(self, "name", "Unnamed Analysis")}'
        plt.title(title)
        plt.grid(True)
        plt.tight_layout()
        plt.show()
        return fig


    def scatter3d(self, dims=(0, 1, 2), annotate=True, n_clusters=None):
        """
        3D scatter plot in selected dimensions with optional cluster-based coloring.

        Parameters
        ----------
        dims : tuple
            Dimensions to plot.
        annotate : bool
            Annotate points with their index.
        n_clusters : int or None
            If provided, use clustering to color points.

        Returns
        -------
        fig : matplotlib.figure.Figure
        """
        fig = plt.figure(figsize=(10, 8))
        ax = fig.add_subplot(111, projection='3d')
        x, y, z = self.coords[:, dims[0]], self.coords[:, dims[1]], self.coords[:, dims[2]]

        if n_clusters is not None:
            labels = self.get_cluster_labels(n_clusters=n_clusters)
            scatter = ax.scatter(x, y, z, c=labels, cmap='tab10', s=50)
        else:
            scatter = ax.scatter(x, y, z, s=50)
        if annotate:
            for i, name in enumerate(self.names):
                ax.text(x[i], y[i], z[i], str(i), fontsize=9)

        ax.set_xlabel(f"Dim {dims[0] + 1}")
        ax.set_ylabel(f"Dim {dims[1] + 1}")
        ax.set_zlabel(f"Dim {dims[2] + 1}")
        title = f'3D PCoA Scatter Plot — {getattr(self, "name", "Unnamed Analysis")}'
        plt.title(title)
        plt.tight_layout()
        plt.show()
        return fig


    def plot_dendrogram(self, truncate_mode=None, p=10):
        """Plot dendrogram from linkage matrix."""
        if self.linkage_matrix is None:
            self.compute_linkage()
        fig=plt.figure(figsize=(10, 6))
        dendrogram(self.linkage_matrix, labels=self.names, truncate_mode=truncate_mode, p=p)
        plt.title("Hierarchical Clustering Dendrogram")
        plt.xlabel("Sample Index or Name")
        plt.ylabel("Distance")
        plt.tight_layout()
        title = f'PCoA Dendrogram — {getattr(self, "name", "Unnamed Analysis")}'
        plt.title(title)
        plt.show()
        return fig

    def cluster(self, t=1.0, criterion='distance'):
        """Assign cluster labels from linkage matrix."""
        if self.linkage_matrix is None:
            self.compute_linkage()
        return fcluster(self.linkage_matrix, t=t, criterion=criterion)

    def best_dimension(self, max_dim=10):
        """Determine optimal dimension by maximizing silhouette score."""
        scores = []
        for d in range(2, min(self.n, max_dim + 1)):
            coords = self.coords[:, :d]
            labels = fcluster(linkage(coords, method='ward'), t=2, criterion='maxclust')
            try:
                score = silhouette_score(coords, labels)
                scores.append((d, score))
            except:
                continue
        if scores:
            best_d, best_score = max(scores, key=lambda x: x[1])
            self.select_dimensions(list(range(best_d)))
            return best_d, best_score
        return None, None

    def save(self, path):
        """Save current analysis to file."""
        with open(path, 'wb') as f:
            pickle.dump(self, f)

    @staticmethod
    def load(path):
        """Load analysis from file."""
        with open(path, 'rb') as f:
            return pickle.load(f)

    def __repr__(self):
        return (f"DNAPairwiseAnalysis(\n"
                f"  name={self.name},\n"
                f"  n_samples={self.n},\n"
                f"  dimensions={len(self.dimensions)},\n"
                f"  active_dims={self.dimensions},\n"
                f"  linkage_computed={'Yes' if self.linkage_matrix is not None else 'No'}\n)")

    def __str__(self):
        return f"DNAPairwiseAnalysis with {self.n} samples in {len(self.dimensions)} dimensions"

    def dimension_variance_curve(self, threshold=0.5, plot=True, figsize=(8, 5)):
        """
        Computes the cumulative explained variance (based on pairwise distances) as a function of
        the number of dimensions used (from 1 to n-1). Optionally plots the curve and the point
        where the threshold (default 0.5) is reached.

        Parameters
        ----------
        threshold : float
            Fraction of total variance to reach (default 0.5).
        plot : bool
            If True, display the variance curve and highlight dhalf.
        figsize : tuple
            Size of the figure if plotted.

        Returns
        -------
        dhalf : int
            Number of dimensions needed to reach the threshold.
        curve : list of float
            Normalized cumulative variance (in [0, 1]) for dimensions 1 to n-1.
        """
        max_d = self.coords.shape[1] # min(self.n - 1, self.coords.shape[1])
        total_var = np.sum(pairwise_distances(self.coords[:, :max_d]))
        curve = []
        for d in range(1, max_d + 1):
            d_var = np.sum(pairwise_distances(self.coords[:, :d]))
            curve.append(d_var / total_var)

        dhalf = next((i + 1 for i, v in enumerate(curve) if v >= threshold), max_d)

        if plot:
            fig=plt.figure(figsize=figsize)
            plt.plot(range(1, max_d + 1), curve, marker='o', label='Cumulative variance')
            plt.axhline(y=threshold, color='gray', linestyle='--', label=f'Threshold = {threshold}')
            plt.axvline(x=dhalf, color='red', linestyle='--', label=f'dhalf = {dhalf}')
            plt.xlabel('Number of dimensions used')
            plt.ylabel('Normalized cumulative variance (distance)')
            title = f'Cumulative Variance vs. Dimensions — {getattr(self, "name", "Unnamed Analysis")}'
            plt.title(title)
            plt.grid(True)
            plt.legend()
            plt.tight_layout()
            plt.show()

        return dhalf, fig


# --------------------------------------------
# DNAsignal_collection class
# This class is intended for 2D intepretation
# --------------------------------------------
class DNAsignal_collection(list):
    """
    A collection of DNAsignal objects (e.g., from a GC-MS chromatogram) supporting symbolic sinusoidal encoding,
    full tensor construction, and blind deconvolution using latent component analysis.

    Purpose
    -------
    DNAsignal_collection is designed to enable symbolic and positional encoding of multiple 1D analytical signals
    (e.g., ion channels from GC-MS data) using lettered segments and sinusoidal encodings. It allows combining
    multiple encoded signals into a single tensor for processing with machine learning methods, including dimensionality
    reduction and blind source separation.

    Theory
    ------
    Each signal is decomposed into symbolic segments based on their local morphology (encoded as letters like A, B, Y, Z, _).
    Each segment is represented in an embedding space of dimension `d` via sinusoidal encoding. The 3D tensor
    $v_{t,m,d}$ is composed of:

    - $E_{t,m,d}$: symbolic encoding across segments.
    - $PE_t$: positional encoding along the time/segment axis (t).
    - $PE_m$: positional encoding along the mass channel/ion axis (m).

    Combining these yields:

    $$
    v_{t,m,d} = E_{t,m,d} + PE_t(t,d) + PE_m(m,d)
    $$

    Methods
    -------
    - `.sinencode_dna_full(scale=4)`: performs symbolic encoding at a given scale.
    - `.E_symbol`: property returning symbolic component E for each letter and scale.
    - `.PE_t`: positional encoding per letter along t for each scale.
    - `.PE_m`: positional encoding along m (mass channels).
    - `.vtm`: dictionary of $v_{t,m,d}$ matrices per letter.
    - `.vtm_full`: complete tensor (sum of all letters) used for machine learning.
    - `.deconvolve_latent_sources(...)`: uses PCA to decompose full tensor into component chromatograms.
    - `.plot_v_symbol_components(...)`: visualizes the construction of $v_{t,m}$ for each letter.
    - `.plot_vtm_full(...)`: visualizes the components of the full $v_{t,m}$ tensor.

    Parameters
    ----------
    signals : list of DNAsignal
        List of DNAsignal instances (e.g., one per ion channel).
    vtmscale : int
        Scale to calculate vtm_full
    rasterscan : bool (default=True)
        True if one point is read at a time by the detector
        In practice, flatten the 2D signal (T × m) into a single 1D time axis by appending
        all temporal channels one after another:
    dtype : type or np.dtype, optional
        Numeric dtype used for storing encoding arrays (E_symbol, PE_t, PE_m, vtm, vtm_full).
        Defaults to np.float32 for reduced memory usage.

    Attributes
    ----------
    m : int
        Number of signals (ion channels).
    d : int
        Embedding dimension.
    letters : list of str
        List of symbolic segment labels used in the encodings.
    scales : list of int
        Available scales for the symbolic encoding.
    _E_symbol : dict
        Cached symbolic encoding tensors for each letter and scale.
    _PE_t : dict
        Cached positional encoding along t for each letter and scale.
    _PE_m : dict
        Cached positional encoding along m for each scale.
    _vtm : dict
        Cached symbolic+positional tensors per letter.
    _vtm_full : np.ndarray
        Cached full encoding tensor combining all letters.
    """

    def __init__(self, *signals, vtmscale=None, rasterscan=True, dtype=np.float32):
        """
        Initialize a DNAsignal_collection from DNAsignal instances.

        Parameters
        ----------
        *signals : DNAsignal
            One or more DNAsignal objects.
        vtmscale : int (default self.scales[0])
            Scale to calculate vtm_full
        dtype : type or np.dtype, optional
            Numeric dtype used for storing encoding arrays (E_symbol, PE_t, PE_m, vtm, vtm_full).
            Defaults to np.float32 for reduced memory usage.
        """
        if not signals:
            raise ValueError("At least one DNAsignal must be provided.")
        if not all(isinstance(sig, DNAsignal) for sig in signals):
            raise TypeError("All elements must be DNAsignal instances.")
        super().__init__(signals)
        self.issinencoded = False
        if not all(hasattr(sig,"codesfull") for sig in signals):
            raise ValueError("All elements must be fully DNA encoded.")

        self.dtype = dtype
        self.d = getattr(signals[0], "d", 128)
        self._PE_t = {}
        self._PE_m = {}
        self._E_symbol = {}
        self._vtm = {}
        self._vtm_full = None
        self.rasterscan = rasterscan

        self.letters = sorted({l for sig in signals for l in getattr(sig, "letters", [])})
        self.scales = sorted({s for sig in signals for s in getattr(sig, "codesfull", {}).keys()})
        if vtmscale is None:
            vtmscale = self.scales[0]
        if not isinstance(vtmscale,int):
            raise TypeError(f"vtmscale must be int not a {type(vtmscale).__name__}")
        if vtmscale not in self.scales:
            raise ValueError(f"vtmscale must be chosen among {[self.scales]}")
        self.vtmscale = vtmscale

    @property
    def m(self):
        """Return the number of signals in the collection."""
        return len(self)

    def sinencode_dna_full(self, d_model=128, N=10000, operation="sum"):
        """
        🌀 Encode all DNAsignal instances using full-resolution sinusoidal embeddings,
        grouped by letter and organized per scale.

        Parameters
        ----------
        d_model : int, optional
            Dimensionality of the sinusoidal embedding (default is 128).
        N : int, optional
            Maximum number of positions for the encoding (default is 10000).
        operation : str or None
            If "sum", sum all position encodings per letter.
            If "mean", average encodings.
            If None, retain full (n_occurrences × d_model) arrays.
        """
        for j in range(self.m):
            self[j].sinencode_dna_full(d_model=d_model, N=N, operation=operation)
        for j in range(self.m):
            self[j].sincodesfull_aggregated.unwrap_letters_to_matrix()
        self.d = d_model
        self.issinencoded = True

    def combine_embeddings(self, selected_letters=None):
        """
        Combine unwrapped embeddings across all signals for each scale.

        Parameters
        ----------
        selected_letters : list of str, optional
            If provided, restrict to these letters.

        Returns
        -------
        dict
            Dictionary {scale: {letter: (m, d)}} for each selected scale and letter.
        """
        results = {}
        for scale in self.scales:
            letters = selected_letters or self.letters
            combined = {letter: np.zeros((self.m, self.d),dtype=self.dtype) for letter in letters}
            for i, sig in enumerate(self):
                if not isinstance(sig,DNAsignal):
                    raise TypeError(f'signal must be a DNAsignal not a {type(sig).__name__}')
                if not hasattr(sig,"sincodesfull_aggregated"):
                    raise ValueError("nothing to combine, call sinencode_dna_full(...) before combining")
                matrix = sig.sincodesfull_aggregated.unwrapped_matrix[scale]
                signal_letters = list(sig.sincodesfull[scale].keys())
                for j, letter in enumerate(signal_letters):
                    if letter in letters:
                        combined[letter][i, :] = matrix[j]
            results[scale] = combined
        return results

    def scale_alignment(self, method="zscore"):
        """
        Normalize embeddings across all signals and all scales.

        Parameters
        ----------
        method : str
            One of {"zscore", "minmax"}.
        """
        for sig in self:
            for scale in self.scales:
                if not isinstance(sig,DNAsignal):
                    raise TypeError(f'signal must be a DNAsignal not a {type(sig).__name__}')
                if not hasattr(sig,"sincodesfull_aggregated"):
                    raise ValueError("nothing to scale, call sinencode_dna_full(...) before scaling")
                mat = sig.sincodesfull_aggregated.unwrapped_matrix[scale]
                if method == "zscore":
                    sig.sincodesfull.unwrapped_matrix[scale] = StandardScaler().fit_transform(mat)
                elif method == "minmax":
                    min_val = mat.min(axis=0)
                    max_val = mat.max(axis=0)
                    sig.sincodesfull.unwrapped_matrix[scale] = (mat - min_val) / (max_val - min_val)
                else:
                    raise ValueError("Unsupported scaling method. Choose 'zscore' or 'minmax'.")

    def reduce_dimensions(self, method="pca", selected_letters=None, n_components=2):
        """
        Apply dimensionality reduction (PCA or UMAP) across signals for each scale.

        Parameters
        ----------
        method : str
            "pca" or "umap".
        selected_letters : list of str, optional
            Restrict to a subset of letters.
        n_components : int
            Number of projection dimensions.

        Returns
        -------
        dict
            Dictionary {scale: ndarray} with shape (m, n_components), one per scale.
        """
        combined_by_scale = self.combine_embeddings(selected_letters)
        projections = {}
        for scale, embeddings in combined_by_scale.items():
            if len(embeddings) == 1:
                # Single letter
                letter = list(embeddings.keys())[0]
                X = embeddings[letter]
            else:
                # Concatenate across letters
                X = np.concatenate([embeddings[letter] for letter in embeddings], axis=1)

            if method == "pca":
                reducer = PCA(n_components=n_components)
            elif method == "umap":
                reducer = UMAP(n_components=n_components)
            else:
                raise ValueError("Method must be 'pca' or 'umap'.")

            projections[scale] = reducer.fit_transform(X)

        return projections

    def __str__(self):
        return f"<DNAsignal_collection: {self.m} signals, {len(self.letters)} letters, {len(self.scales)} scales, d={self.d}>"

    def __repr__(self):
        """
        Return a readable summary of the DNAsignal_collection contents.
        Shows the number of signals, available letters, scales, and embedding dimension.
        Flags whether sinusoidal encoding has been performed.
        """
        flag = "✔️" if self.issinencoded else "⚠️ not encoded"
        print(f"<DNAsignal_collection: {self.m} signals>",
              f"\tletters={self.letters}",
              f"\tscales={self.scales}",
              f"\td={self.d}",
              f"\t{'sin-encoded' if flag else 'sin-encoding pending'}",sep="\n"
            )
        return self.__str__()

    def to_dataframe(self, selected_letters=None):
        """
        Export combined embeddings for all scales as a tidy pandas DataFrame,
        suitable for machine learning tasks.

        Parameters
        ----------
        selected_letters : list of str, optional
            Subset of letters to include. If None, include all letters.

        Returns
        -------
        pd.DataFrame
            A long-form DataFrame with columns:
            ['signal_index', 'scale', 'letter', 'dim_0', ..., 'dim_{d-1}']
        """
        import pandas as pd

        all_records = []
        combined_all = self.combine_embeddings(selected_letters)

        for scale, embeddings in combined_all.items():
            for letter, mat in embeddings.items():
                for i, vec in enumerate(mat):
                    record = {"signal_index": i, "scale": scale, "letter": letter}
                    record.update({f"dim_{j}": v for j, v in enumerate(vec)})
                    all_records.append(record)

        return pd.DataFrame(all_records)

    # Generic plot
    def plot(self, letters = None, scales = None, figsize=(18, 10), max_legend=25):
        """
        Plot the encoded signals in subplots.
        Rows represent letters, columns represent scales.
        Each subplot contains overlaid colored curves from all signals.

        Parameters
        ----------
        letters : list or None
            Letters to be plotted. If None, all available letters are plotted.
        scales : list or None
            Scales to be plotted. If None, all available scales are plotted.
        figsize : tuple
            Size of the full figure.
        max_legend : int
            Maximum number of signals to label in the legend.

        Returns
        -------
        matplotlib.figure.Figure
            The figure containing the plots.
        """
        if not self.issinencoded:
            raise ValueError("Signals must be sinencoded first (call sinencode_dna_full).")

        letters = letters or sorted(self.letters)
        scales = scales or sorted(self.scales)
        m = self.m
        d = self.d

        fig, axs = plt.subplots(len(letters), len(scales), figsize=figsize, squeeze=False)
        cmap = cm.get_cmap("plasma", m)

        for i, letter in enumerate(letters):
            for j, scale in enumerate(scales):
                ax = axs[i, j]
                for idx, sig in enumerate(self):
                    if not isinstance(sig,DNAsignal):
                        raise TypeError(f'signal for letter "{letter}" must be a DNAsignal not a {type(sig).__name__}')
                    if not hasattr(sig,"sincodesfull_aggregated"):
                        raise ValueError("nothing to plot, call sinencode_dna_full(...) before plotting")
                    try:
                        vec = sig.sincodesfull_aggregated[scale][letter]
                        if isinstance(vec, np.ndarray) and vec.ndim == 2:
                            vec = vec.mean(axis=0)
                        ax.plot(vec, label=f"sig {idx}", color=cmap(idx / max(1, m-1)))
                    except KeyError:
                        continue
                ax.set_title(f"Letter '{letter}' — Scale {scale}", fontsize=10)
                ax.set_xticks([])
                ax.set_yticks([])
                if i == len(letters) - 1:
                    ax.set_xlabel("Embedding dim")
                if j == 0:
                    ax.set_ylabel("Amplitude")

        # Legend outside of the plot
        handles, labels = axs[0, 0].get_legend_handles_labels()
        if handles:
            fig.legend(handles[:max_legend], labels[:max_legend],
                       loc='upper center', bbox_to_anchor=(0.5, 1.03),
                       ncol=min(max_legend, 5), fontsize="small", title="Signal Index")

        plt.tight_layout()
        return fig

    def plot_letters(self, scale=None, figsize=(12, 6), cmap="viridis"):
        """
        Plot a heatmap of the letter codes (symbolic DNA) across all signals.

        Parameters
        ----------
        scale : int, optional
            Scale to use (default: self.vtmscale).
        figsize : tuple
            Size of the figure.
        cmap : str
            Matplotlib colormap name (default: "viridis").

        Returns
        -------
        matplotlib.figure.Figure
        """
        scale = scale or self.vtmscale
        m = self.m
        T = len(self[0].codesfull[scale])
        letters = self.letters
        letter_to_int = {l: i for i, l in enumerate(letters)}
        int_to_letter = {i: l for l, i in letter_to_int.items()}

        # Build numeric matrix
        mat = np.zeros((m, T), dtype=int)
        counts = {l: 0 for l in letters}
        for i, sig in enumerate(self):
            seq = sig.codesfull[scale]
            for t, c in enumerate(seq):
                mat[i, t] = letter_to_int.get(c, 0)
                counts[c] = counts.get(c, 0) + 1

        # Plot
        fig, ax = plt.subplots(figsize=figsize)
        im = ax.imshow(mat, aspect='auto', cmap=plt.get_cmap(cmap, len(letters)))

        # Legend
        patches = [mpatches.Patch(color=im.cmap(i), label=f"{l} ({counts[l]})") for i, l in enumerate(letters)]
        ax.legend(handles=patches, bbox_to_anchor=(1.05, 1), loc='upper left', title="Letters")

        ax.set_xlabel("Time index")
        ax.set_ylabel("Signal index")
        ax.set_title(f"Letter Code Heatmap (T×m = {T}×{m}) – Occurrences: " +
                     ", ".join(f"{l}:{counts[l]}" for l in letters))

        plt.tight_layout()
        return fig


    def plot_embedding_projection(self, letters=None, scales=None, method='pca', max_points=25, figsize=(14, 10)):
        """
        Plot embedding projections of the encoded signals using PCA (default) or other DR methods.

        Parameters
        ----------
        collection : DNAsignal_collection
            The collection of encoded signals.
        scales : list or None
            Scales to include in the projection. If None, all available scales are used.
        method : str
            Dimensionality reduction method ('pca' only supported for now).
        max_points : int
            Maximum number of signal points to label explicitly.
        figsize : tuple
            Size of the figure.

        Returns
        -------
        fig : matplotlib.figure.Figure
        """
        # Determine scales
        if scales is None:
            scales = self.scales

        letters = letters or sorted(self.letters)
        n_rows = len(letters) + 1  # +1 for the full embedding
        n_cols = len(scales)

        fig, axs = plt.subplots(n_rows, n_cols, figsize=(figsize[0], figsize[1]))
        axs = np.array(axs).reshape((n_rows, n_cols))

        norm = mcolors.Normalize(vmin=0, vmax=self.m)
        cmap = cm.get_cmap("plasma", self.m)

        for j, scale in enumerate(scales):
            combined = self.combine_embeddings()[scale]
            for i, letter in enumerate(letters):
                embeddings = combined[letter]
                reducer = PCA(n_components=2)
                proj = reducer.fit_transform(embeddings)

                ax = axs[i, j]
                for idx, point in enumerate(proj):
                    color = cmap(norm(idx))
                    ax.scatter(*point, color=color)
                    if idx < max_points:
                        ax.text(point[0], point[1], str(idx), fontsize=8)
                ax.set_title(f"Letter '{letter}', Scale {scale}")
                ax.set_xticks([])
                ax.set_yticks([])
                ax.grid(True)

            # Last row: combined embeddings for all letters
            all_embeddings = [combined[letter] for letter in letters]
            all_embeddings = np.concatenate(all_embeddings, axis=1)
            reducer = PCA(n_components=2)
            proj = reducer.fit_transform(all_embeddings)
            ax = axs[-1, j]
            for idx, point in enumerate(proj):
                color = cmap(norm(idx))
                ax.scatter(*point, color=color)
                if idx < max_points:
                    ax.text(point[0], point[1], str(idx), fontsize=8)
            ax.set_title(f"All Letters, Scale {scale}")
            ax.set_xticks([])
            ax.set_yticks([])
            ax.grid(True)

        fig.suptitle("Embedding Projections by Letter and Scale", fontsize=16)
        plt.tight_layout()
        return fig

    # @property: def E_symbol(self) -> Dict[int, Dict[str, np.ndarray]]:
    @property
    def E_symbol(self):
        """
        E_symbol(t, m): Symbolic component of the 2D encoding for each scale and letter.

        For each letter, builds a tensor (T_letter, M, D) where:
        - T_letter is the total number of segments of that letter across all M signals,
        - D is the embedding dimension,
        - M is the number of ion channels.

        Returns
        -------
        dict of dict
            Mapping {scale: {letter: ndarray of shape (T_letter, m, d)}}
        """
        if not self.issinencoded:
            raise ValueError("Call sinencode_dna_full() first.")
        if not all(isinstance(sig, DNAsignal) for sig in self):
            raise TypeError("E_symbol is only defined for DNAsignal_collection containing DNAsignal instances.")
        if not all(hasattr(sig, "sincodesfull") for sig in self):
            raise ValueError("Missing `sincodesfull` on at least one DNAsignal. Run `.sinencode_dna_full()` first.")
        if self._E_symbol:
            return self._E_symbol

        self._E_symbol = {}
        for scale in self.scales:
            E_scale = {}
            for letter in self.letters:
                # First pass: compute total segments
                segment_counts = [
                    sig.sincodesfull[scale][letter].shape[0]
                    if letter in sig.sincodesfull[scale] else 0
                    for sig in self
                ]
                total_segments = sum(segment_counts)

                if total_segments == 0:
                    continue  # Skip completely missing letter at this scale

                # Allocate the global array
                E_mat = np.zeros((total_segments, self.m, self.d),dtype=self.dtype)
                offset = 0
                for i, sig in enumerate(self):
                    segments = sig.sincodesfull[scale].get(letter)
                    if segments is None or segments.size == 0:
                        continue
                    n_segments = segments.shape[0]
                    E_mat[offset:offset + n_segments, i, :] = segments
                    offset += n_segments
                E_scale[letter] = E_mat
            self._E_symbol[scale] = E_scale
        return self._E_symbol

    # @property: def PE_t(self) -> Dict[int, Dict[str, np.ndarray]]:
    @property
    def PE_t(self):
        """
        PE_t(t): Positional encoding along t (segment axis) per scale and letter.

        For each scale and letter, this provides a matrix of shape (n_segments, d),
        where n_segments is the total number of segments of that letter across the m signals,
        and d is the embedding dimension.

        Returns
        -------
        dict of dict
            Mapping {scale: {letter: ndarray of shape (n_segments, d)}}
        """
        if not self.issinencoded:
            raise ValueError("Call sinencode_dna_full() first.")
        if not all(isinstance(sig, DNAsignal) for sig in self):
            raise TypeError("PE_t is only defined for DNAsignal_collection containing DNAsignal instances.")
        if not all(hasattr(sig, "sincodesfull") for sig in self):
            raise ValueError("Missing `sincodesfull` on at least one DNAsignal. Run `.sinencode_dna_full()` first.")
        if self._PE_t:
            return self._PE_t

        self._PE_t = {}
        for scale in self.scales:
            PE_t_scale = {}
            for letter in self.letters:
                # First pass: compute total segments
                total_segments = sum(
                    sig.sincodesfull[scale][letter].shape[0]
                    if letter in sig.sincodesfull[scale] else 0
                    for sig in self
                )

                if total_segments == 0:
                    continue  # Skip completely missing letter at this scale

                # Generate sinusoidal encoding
                position = np.arange(total_segments)[:, np.newaxis]  # (T_letter, 1)
                div_term = np.linspace(0, 1, self.d)[np.newaxis, :]  # (1, d)
                PE_t_scale[letter] = np.sin(position * div_term)     # (T_letter, d)

            self._PE_t[scale] = PE_t_scale
        return self._PE_t


    # @property: def PE_m(self) -> Dict[int, np.ndarray]:
    @property
    def PE_m(self):
        """
        PE_m(m): Positional encoding along m (IC index axis) per scale.

        Returns
        -------
        dict
            Mapping: {scale: array of shape (m, d)}
        """
        if not self.issinencoded:
            raise ValueError("Call sinencode_dna_full() first.")
        if not all(isinstance(sig, DNAsignal) for sig in self):
            raise TypeError("PE_m is only defined for DNAsignal_collection containing DNAsignal instances.")
        if not all(hasattr(sig, "sincodesfull") for sig in self):
            raise ValueError("Missing `sincodesfull` on at least one DNAsignal. Run `.sinencode_dna_full()` first.")
        if self._PE_m:
            return self._PE_m
        for scale in self.scales:
            PE_m_scale = np.zeros((self.m, self.d),dtype=self.dtype)
            count = 0
            for i in range(self.m):
                for letter in self.letters:
                    if letter in self[i].sincodesfull[scale]:
                        PE_m_scale[i] += self[i].sincodesfull_aggregated[scale][letter]
                        count += 1
            self._PE_m[scale] = PE_m_scale / len(self.letters)
        return self._PE_m

    @property
    def vtm(self):
        """
        Compute the full encoded matrix v_{t,m} for each letter at each scale.

        This property combines three orthogonal components:
        - E_symbol(t, m): the original per-segment encoding for each letter
        - PE_t(t): a sinusoidal encoding applied along the segment (time) axis
        - PE_m(m): a sinusoidal encoding applied along the signal (IC) axis

        The resulting tensor for each scale and letter is of shape (n_segments, m, d), where:
            - n_segments: number of segments (time positions) per letter
            - m: number of DNAsignal instances in the collection
            - d: dimensionality of the encoding space (d_model)

        Returns
        -------
        dict
            A nested dictionary {scale: {letter: array of shape (n_segments, m, d)}}.
        """
        if self._vtm:
            return self._vtm
        if not self.issinencoded:
            raise ValueError("Call sinencode_dna_full() first.")
        vtm_result = {}
        for scale in self.scales:
            result = {}
            PE_m_seg = self.PE_m[scale]
            for letter in self.letters:
                E_seg = self.E_symbol[scale][letter]
                PE_t_seg = self.PE_t[scale][letter]
                PE_t_broadcast = PE_t_seg[:, np.newaxis, :]  # (n_segments, 1, d)
                PE_m_broadcast = PE_m_seg[np.newaxis, :, :]  # (1, m, d)

                if E_seg.shape[0] != PE_t_broadcast.shape[0]:
                    raise ValueError(f"Mismatch: E_symbol[{scale}][{letter}].shape[0]={E_seg.shape[0]} "
                                     f"but PE_t[{scale}][{letter}].shape[0]={PE_t_broadcast.shape[0]}")
                v_tm = E_seg + PE_t_broadcast + PE_m_broadcast
                result[letter] = v_tm
            vtm_result[scale] = result
        self._vtm = vtm_result
        return self._vtm

    def plot_v_symbol_components(self, scale=None, dims="all"):
        """
        Plot the components E_symbol, PE_t, PE_m and their sum v_{t,m}
        as image matrices for each letter at a given scale.

        Parameters
        ----------
        scale : int, optional
            The scale to use. Defaults to the first available scale.
        dims : "all", list or slice
            Which dimensions to include in the sum over `d`. Default is all.

        Returns
        -------
        matplotlib.figure.Figure
        """
        if not self.issinencoded:
            raise ValueError("Call sinencode_dna_full() first.")

        if scale is None:
            scale = self.scales[0]
        if scale not in self.scales:
            raise ValueError(f"Scale {scale} not found in collection.")

        fig, axes = plt.subplots(len(self.letters), 4, figsize=(16, 4 * len(self.letters)))
        fig.suptitle(f"Scale {scale}: Components of $v_{{t,m}} = E + PE_t + PE_m$", fontsize=14)

        PE_m = self.PE_m[scale]
        for i, letter in enumerate(self.letters):
            E = self.E_symbol[scale][letter]
            PE_t = self.PE_t[scale][letter][:, np.newaxis, :]
            PE_m_broad = PE_m[np.newaxis, :, :]
            V = self.vtm[scale][letter]

            # Validate dims
            d = E.shape[-1]
            if dims == "all":
                dims_to_use = np.arange(d)
            elif isinstance(dims, (list, slice, np.ndarray)):
                dims_to_use = np.arange(d)[dims]
            else:
                raise TypeError("`dims` must be 'all', a list, or a slice.")

            # Subset and sum
            E_2D = E[:, :, dims_to_use].sum(axis=-1)
            PE_t_2D = PE_t[:, :, dims_to_use].sum(axis=-1)
            PE_m_2D = PE_m_broad[:, :, dims_to_use].sum(axis=-1)
            V_2D = V[:, :, dims_to_use].sum(axis=-1)

            for ax, mat, title in zip(axes[i], [E_2D, PE_t_2D, PE_m_2D, V_2D],
                                      [f"E_symbol[{letter}]", f"PE_t[{letter}]", f"PE_m[{letter}]", f"v_t,m[{letter}]"]):
                im = ax.imshow(mat, aspect='auto', cmap='viridis')
                ax.set_title(title)
                fig.colorbar(im, ax=ax)

        return fig

    # -- full interpretation for 2D signals ---
    @property
    def vtm_full(self):
        """
        Compute and store the full encoded tensor for the GC-MS signal:
            - If self.rasterscan is False:
                shape = (T, m, d)
            - If self.rasterscan is True:
                shape = (T * m, d)

        This combines:
        - Symbolic embedding per character (one-hot or learned)
        - Positional encoding along time axis
        - PE_m (mass/IC identity) is used only if rasterscan is False

        Returns
        -------
        np.ndarray
            Encoded array (T, m, d) or (T*m, d)
        """
        if hasattr(self, "_vtm_full") and self._vtm_full is not None:
            return self._vtm_full

        if not self.issinencoded:
            raise ValueError("Call .sinencode_dna_full() first.")

        if not all(isinstance(sig, DNAsignal) for sig in self):
            raise TypeError("vtm_full is only defined for DNAsignal_collection containing DNAsignal instances.")

        scale = self.vtmscale
        m = self.m
        d = self.d
        T = len(self[0].codesfull[scale])  # All signals assumed same length

        # -- Build symbolic embedding per letter (1-hot or learned) --
        letters = self.letters
        letter_indices = {c: i for i, c in enumerate(letters)}
        rng = np.random.default_rng(42)
        embedding_matrix = rng.standard_normal((len(letters), d)).astype(self.dtype)
        embedding_dict = {c: embedding_matrix[i] for i, c in enumerate(letters)}

        # -- Build encoded E_list --
        E_list = []
        for sig in self:
            code_str = sig.codesfull[scale]
            segments = np.stack([embedding_dict[c] for c in code_str], axis=0)  # (T, d)
            E_list.append(segments)

        # -- Positional encoding along time --
        PE_t = SinusoidalEncoder(d_model=d, N=T).encode(np.arange(T))  # (T, d)

        if self.rasterscan:
            # -- Flatten all signals: shape = (T * m, d)
            E_flat = np.concatenate(E_list, axis=0)  # (T * m, d)
            PE_t_flat = np.tile(PE_t, (m, 1))  # repeat time encoding m times
            self._vtm_full = E_flat + PE_t_flat
        else:
            # -- Stack signals into (T, m, d)
            E_tensor = np.stack(E_list, axis=1)  # (T, m, d)
            PE_t_tensor = PE_t[:, np.newaxis, :]  # (T, 1, d)
            PE_m_tensor = self.PE_m[scale][np.newaxis, :, :]  # (1, m, d)
            self._vtm_full = E_tensor + PE_t_tensor + PE_m_tensor  # full positional tensor

        return self._vtm_full


    def plot_vtm_full(self, scale=None, dims="all"):
        """
        Visualize the components and sum of the full encoded GC-MS signal at a given scale.

        Parameters
        ----------
        scale : int, optional
            Scale to visualize. Defaults to the first scale.
        dims : "all", list or slice
            Dimensions of the embedding `d` to include (default: all).

        Returns
        -------
        matplotlib.figure.Figure
        """
        if not self.issinencoded:
            raise ValueError("Call sinencode_dna_full() first.")

        scale = scale or self.vtmscale
        if scale not in self.scales:
            raise ValueError(f"Scale {scale} is not valid. Available scales: {self.scales}")

        v_tm = self.vtm_full
        if self.rasterscan:
            Tm, d = v_tm.shape
            T = len(self[0].codesfull[scale])
            m = len(self)
            v_tm = v_tm.reshape(T, m, d)

        T, m, d = v_tm.shape

        if dims == "all":
            dims_to_use = np.arange(d)
        elif isinstance(dims, (list, slice, np.ndarray)):
            dims_to_use = np.arange(d)[dims]
        else:
            raise TypeError("`dims` must be 'all', a list, or a slice.")

        PE_t = SinusoidalEncoder(d_model=d, N=T).encode(np.arange(T))[:, dims_to_use]
        E_partial = v_tm[:, :, dims_to_use] - PE_t[:, None, :]

        if self.rasterscan:
            PE_m_2D = np.zeros_like(E_partial.sum(axis=2))  # dummy for display
        else:
            PE_m = self.PE_m[scale][:, dims_to_use]
            E_partial -= PE_m[None, :, :]
            PE_m_2D = PE_m.sum(axis=1)[None, :].repeat(T, axis=0)

        E_2D = E_partial.sum(axis=2)
        PE_t_2D = PE_t.sum(axis=1)[:, None].repeat(m, axis=1)
        V_2D = v_tm[:, :, dims_to_use].sum(axis=2)

        fig, axs = plt.subplots(1, 4, figsize=(20, 5))
        titles = ["E_symbol total", "PE_t", "PE_m" if not self.rasterscan else "PE_m (omitted)", "v_{t,m} total"]
        matrices = [E_2D, PE_t_2D, PE_m_2D, V_2D]

        for ax, mat, title in zip(axs, matrices, titles):
            im = ax.imshow(mat, aspect='auto', cmap='viridis')
            ax.set_title(title)
            ax.set_xlabel("m")
            ax.set_ylabel("t")
            fig.colorbar(im, ax=ax)

        fig.suptitle(f"Full GC-MS Representation (sum over dims={dims}) at scale={scale}", fontsize=14)
        plt.tight_layout()
        return fig


    def deconvolve_latent_sources(self, n_components=64, inertia_loss_threshold=0.25, plot=True, nmax_plot=8):
        """
        Perform dimensionality reduction on the 3D tensor v_{t,m,d} to extract non-coeluted
        compound chromatograms using PCA, with optional plotting.

        Parameters
        ----------
        n_components : int, optional
            Maximum number of latent components (e.g., pure compounds) to extract.
            Default is 64.
        inertia_loss_threshold : float, optional
            The maximum allowed proportion of total variance to lose in the projection.
            Default is 0.25 (i.e., at least 75% of the variance should be preserved).
        plot : bool, optional
            Whether to display diagnostic plots.
        nmax_plot : int, optional
            Maximum number of components to visualize in plots.

        Returns
        -------
        components : np.ndarray
            Component matrix of shape (n_selected_components, D), representing
            the spectral basis vectors (latent features).
        chromatograms : np.ndarray
            Projected chromatograms for each component, shape (T, M, n_selected_components).
        explained_variance_ratio : np.ndarray
            Variance explained by each selected component.
        """

        def find_pca_corner(cumulative_variance, min_index=0):
            """
            Identify the 'corner' in a monotonic increasing curve using the max
            deviation from the straight line between cumulative_variance[min_index]
            and the last point.
            Returns: corner_index or None
            """
            x = np.arange(len(cumulative_variance))
            x0, x1 = x[min_index], x[-1]
            y0, y1 = cumulative_variance[min_index], cumulative_variance[-1]
            # Linear baseline from (x0, y0) to (x1, y1)
            baseline = y0 + (x - x0) * (y1 - y0) / (x1 - x0)
            # Deviation from baseline
            delta = cumulative_variance - baseline
            corner_index = np.argmax(delta)-1 # -1 since 0 has been added
            # Diagnostic check: is the corner meaningful?
            post_corner_trend = np.mean(np.diff(cumulative_variance[corner_index:]))
            if corner_index>0 and delta[corner_index] < 2 * post_corner_trend:
                return None  # not significant enough
            return corner_index


        vtm_full = self.vtm_full
        if self.rasterscan:
            T = len(self[0].codesfull[self.vtmscale])
            m = len(self)
            Tm, D = vtm_full.shape
            vtm_full = vtm_full.reshape(T, m, D)

        T, M, D = vtm_full.shape
        X = vtm_full.reshape(-1, D)

        pca = PCA(n_components=min(n_components, D))
        X_transformed = pca.fit_transform(X)
        explained_variance_ratio = pca.explained_variance_ratio_
        cumulative_variance = np.insert(np.cumsum(explained_variance_ratio), 0, 0.0)

        # select PCA components based on inertia_loss_threshold
        n_selected = np.searchsorted(cumulative_variance, 1 - inertia_loss_threshold)
        # refinement based on corners
        n_corner = find_pca_corner(cumulative_variance)
        n_selected = n_corner if n_corner is not None else n_selected

        # optimization for plotting
        n_selected = min(n_selected, n_components) # requested
        n_plot = min(nmax_plot, n_selected) # required for plotting
        nrow = math.ceil(math.sqrt(n_plot))
        ncol = math.ceil(n_plot / nrow)  # needed for filling the tiles
        n_plot = min(n_components,max(n_selected,nrow*ncol))
        n_selected = max(n_selected,n_plot) if plot else n_selected

        selected_scores = X_transformed[:, :n_selected]
        chromatograms = selected_scores.reshape(T, M, n_selected)
        components = pca.components_[:n_selected]

        if plot:
            fig1, ax = plt.subplots(figsize=(8, 4))
            ax.plot(np.arange(0, len(cumulative_variance)), cumulative_variance, marker='o')
            if n_corner:
                ax.axvline(n_corner, color='DeepSkyBlue', linestyle='--', label=f'Corner: {n_corner}')
            ax.axvline(n_selected, color='Tomato', linestyle='--', label=f'Selected: {n_selected}')
            ax.set_title("Cumulative Explained Variance")
            ax.set_xlabel("Number of Components")
            ax.set_ylabel("Cumulative Variance")
            ax.legend()
            plt.tight_layout()
            plt.show()

            nrow = math.ceil(math.sqrt(n_plot)) # updated
            ncol = math.ceil(n_plot / nrow)
            fig2, axs = plt.subplots(nrow, ncol, figsize=(4 * ncol, 4 * nrow))
            axs = axs.flatten()
            for i in range(min(n_plot,nrow*ncol)):
                data = chromatograms[:, :, i]
                im = axs[i].imshow(data, aspect='auto', cmap='plasma')
                axs[i].set_title(f"[{i + 1}/{n_components}] 2D signal: ({explained_variance_ratio[i]*100:.3g} %)")
                axs[i].set_xlabel("dim2: e.g. Mass (m or m/z)")
                axs[i].set_ylabel("dim1: e.g. Time (t)")
                fig2.colorbar(im, ax=axs[i])
            plt.tight_layout()
            plt.show()

            fig3, axs = plt.subplots(n_plot, 2, figsize=(10, 3 * n_plot))
            for i in range(n_plot):
                data = chromatograms[:, :, i]
                axs[i][0].plot(data.sum(axis=1))
                axs[i][0].invert_yaxis()
                axs[i][0].set_title(f"[{i + 1}/{n_components}] dim1 projection: e.g. TIC ({explained_variance_ratio[i]*100:.3g} %)")
                axs[i][1].plot(data.sum(axis=0))
                axs[i][1].set_title(f"[{i + 1}/{n_components}] dim2 projection: e.g. m/z spectrum ({explained_variance_ratio[i]*100:.3g} %)")
            plt.tight_layout()
            plt.show()

            fig4, axs = plt.subplots(nrows=n_plot, figsize=(12, 2.5 * n_plot), sharex=True)
            if n_plot == 1:
                axs = [axs]
            for i in range(n_plot):
                axs[i].plot(components[i])
                axs[i].set_title(f"PCA Component #{i + 1} ({explained_variance_ratio[i]*100:.3g} %)")
            plt.tight_layout()
            plt.show()

            figs = {"variance":fig1,"chromatograms":fig2,
                    "projections":fig3,"PCAcomponents":fig4}

        else:
            figs = None

        return components, chromatograms, explained_variance_ratio[:n_selected], figs

# %% Signal and signal_collection classes

# ------------------------
# signal
# ------------------------
class signal:
    """
    signal: A self-documented 1D analytical signal container for reproducible scientific workflows.

    This class is designed for lab-grade signal processing and traceable data storage.
    It represents a discrete 1D signal (e.g., chromatogram, spectrum, transient) with full metadata,
    support for symbolic transformation, numerical operations, plotting, and structured saving/loading.

    Key features include:
    - Portable metadata (user, time, host, cwd, version)
    - Domain-aware plots and operations
    - Reproducible signal serialization in JSON or compressed format
    - Full traceability of all transformation events
    - Optional recursive backup of prior states

    Attributes
    ----------
    x : np.ndarray
        Sampling domain (e.g., time, wavelength, chemical shift).
    y : np.ndarray
        Signal values aligned with x.
    name : str
        Label for plots and file storage (used as default filename).
    type : str
        Optional tag (e.g., 'GC-MS', 'FTIR', 'NMR', 'synthetic').
    x_label : str
        Label for the x-axis (e.g., 'wavenumber').
    x_unit : str
        Unit of the x-axis (e.g., 'cm⁻¹').
    source : str
        Origin label ('array', 'peaks', 'noise', 'imported'...).
    metadata : dict
        Includes user, date, host, cwd, version — filled automatically unless overridden.
    color (str or [rgb]), linestyle (str), linewidth (str)
    _previous: signal
        deep-copy of current object
    _history: dict
        "user@host:timestamp | uidkey" :{"action":str, "details": str}

    Key Methods
    -----------
    - normalize(...)          : Normalize the signal to positive values
    - from_peaks(...)         : Construct signal from a `peaks` object
    - add_noise(...)          : Return noisy variant (Poisson, Gaussian, ramp or constant bias)
    - align_with(...)         : Align this signal with another (same x domain)
    - copy()                  : Deep copy
    - save(...)               : Save as JSON or .gz (optional CSV export)
    - load(...)               : Load from saved file
    - plot(...)               : Plot the signal with axis labels
    - backup(...)             : Backup current signal (deep-copy stored in _previous)
    - restore(...)            : Restore the previous state of the signal
    - apply_poisson_baseline_filter(...) : Apply a Poisson-based filter
    - enable_fullhistory      : enable full history
    - disable_fullhistory     : disable full history
    - _toDNA(signal)          : DNAsignal

    Overloaded Operators
    --------------------
    - +, -, *, /              : Operates on signals or scalars, aligns if needed
    - +=, -=, *=, /=          : In-place functional versions (returns new signal)

    Low-level Methods
    -----------------
    - _current_stamp()        : stamp for events (static method)
    - _copystatic()           : deep-copy of signal only (use copy for a full copy) (static method)
    - _events()               : register a processing step
    - _to_serializable        : Convert the signal into a dictionary suitable for JSON export
    - _from_serizalizable     : convert a dict (e.g., from JSON import) to signal

    Example
    -------
    >>> s = signal(x, y, name="sample", type="FTIR", x_label="wavenumber", x_unit="cm⁻¹")
    >>> s.add_noise("gaussian", 0.05).plot()
    >>> s.save()  # saves to ./sample.json.gz
    >>> s2 = signal.load("sample.json.gz")
    """

    def __init__(self, x=None, y=None,
                 name="signal",
                 type="generic",
                 x_label="index",x_unit="-",
                 y_label="intensity",y_unit="a.u.",
                 metadata=None,
                 source="array",
                 user=None, date=None, host=None, cwd=None, version=None,
                 color=None,linewidth=2,linestyle='-',message=None,fullhistory=True,):
        """
        Initialize a signal instance.

        Parameters
        ----------
        x : array-like, optional
            Sampling domain (e.g., time, wavelength, chemical shift).
        y : array-like, optional
            Signal values aligned with x.
        name : str, optional
            Label for plots and file storage (used as default filename).
        type : str, optional
            Optional tag (e.g., 'GC-MS', 'FTIR', 'NMR', 'synthetic').
        x_label : str, optional
            Label for the x-axis (e.g., 'wavenumber').
        x_unit : str, optional
            Unit of the x-axis (e.g., 'cm⁻¹').
        y_label : str, optional
            Label for the y-axis (e.g., 'intensity').
        y_unit : str, optional
            Unit of the y-axis (e.g., 'a.u.').
        source : str, optional
            Origin label ('array', 'peaks', 'noise', 'imported'...).
        metadata : dict, optional
            Full metadata dictionary. If None, fields below are auto-filled.
        user : str, optional
            Username (defaults to current user).
        date : str, optional
            ISO timestamp of creation.
        host : str, optional
            Hostname of the machine.
        cwd : str, optional
            Current working directory at creation time.
        version : str, optional
            Software version (defaults to global __version__).
        color (default=None), linestyle (default="-"), linewidth (default=2): optional
            Plot styling.
        message : str, optional
            message to record
        fullhistory : bool, optional (default=True)
            flag to enable full history
        """

        self.x = np.asarray(x) if x is not None else None
        self.y = np.asarray(y) if y is not None else None
        self.name = name
        self.type = type
        self.x_label = x_label
        self.x_unit = x_unit
        self.y_label = y_label
        self.y_unit = y_unit
        self.source = source
        self.metadata = metadata or {
            "user": user or getpass.getuser(),
            "date": date or datetime.datetime.now().isoformat(),
            "host": host or socket.gethostname(),
            "cwd": cwd or os.getcwd(),
            "version": version or globals().get("__version__", "undefined"),
            "other": ""
        }
        self.color = color
        self.linestyle = linestyle
        self.linewidth = linewidth
        self._previous = None
        self._history = {}
        self._fullhistory = fullhistory  # Do not serialize this field
        self._events("init", {"from": self.source, "message": message})

    def normalize(self, mode="zscore+shift", inplace=True, shift_eps=1e-6):
        """
        Normalize the signal to positive values using different normalization strategies.

        Parameters
        ----------
        mode : str
            Normalization mode:
            - "zscore+shift"   : (y - mean) / std, then shift so min is shift_eps
            - "minmax"         : (y - min) / (max - min), scales to [0, 1]
            - "max"            : y / max, scales to [0, 1]
            - "l1"             : y / sum(|y|), sums to 1 (like probability)
            - "energy"         : y / sqrt(sum(y^2)), unit energy
            - "none"           : No normalization, just returns a copy or itself

        inplace : bool
            Whether to modify the signal in place. If False, returns a new signal.

        shift_eps : float
            Minimum value to add after z-score shift to ensure strictly positive output.

        Returns
        -------
        signal or None
            Normalized signal (if inplace is False), else None.

        Raises
        ------
        ValueError
            If the normalization fails (e.g., due to division by zero).
        """
        if self.y is None:
            raise ValueError("Signal values (y) are not set.")

        y = self.y.astype(float)

        if mode == "zscore+shift":
            mean, std = np.mean(y), np.std(y)
            if std == 0:
                raise ValueError("Zero std deviation: cannot z-score normalize.")
            y = (y - mean) / std
            y -= np.min(y)  # ensure min is zero
            y += shift_eps  # ensure strictly positive
        elif mode == "minmax":
            miny, maxy = np.min(y), np.max(y)
            if maxy - miny == 0:
                raise ValueError("Constant signal: cannot min-max normalize.")
            y = (y - miny) / (maxy - miny)
        elif mode == "max":
            maxy = np.max(np.abs(y))
            if maxy == 0:
                raise ValueError("Signal is zero everywhere.")
            y = y / maxy
        elif mode == "l1":
            total = np.sum(np.abs(y))
            if total == 0:
                raise ValueError("L1 norm is zero.")
            y = y / total
        elif mode == "energy":
            energy = np.sqrt(np.sum(y ** 2))
            if energy == 0:
                raise ValueError("Zero energy: cannot normalize.")
            y = y / energy
        elif mode == "none":
            y = y.copy()
        else:
            raise ValueError(f"Unknown normalization mode: {mode}")

        if inplace:
            self.y = y
            self._events("normalize", {"mode": mode})
        else:
            s = self.copy()
            s.y = y
            s._events("normalize", {"mode": mode})
            return s

    @property
    def n(self):
        """Return the length of the signal and None if it is None"""
        return len(self.y) if self.y is not None else None

    @staticmethod
    def _current_stamp():
        timestamp = datetime.datetime.now().isoformat(timespec="seconds")
        userhost = f"{getpass.getuser()}@{socket.gethostname()}"
        idstamp = uuid.uuid4().hex[:6]
        return f"{timestamp}:{userhost} | {idstamp}"

    @staticmethod
    def _copystatic(s):
        """Static copy of a signal (x and y only)"""
        return s.__class__(x=s.x.copy(), y=s.y.copy(), name=s.name + "_backup")

    def _events(self, action, details=None):
        """
        Register a traceable action in the signal's history.

        Each event is stored with a unique timestamp-based key and includes:
        - user@host
        - timestamp
        - action (string)
        - details (optional dictionary)

        Parameters
        ----------
        action : str
            Description of the event (e.g., 'baseline_filter', 'restore').
        details : dict, optional
            Additional parameters relevant to the action.
        """
        if not hasattr(self, "_history") or not isinstance(self._history, dict):
            self._history = {}
        key = signal._current_stamp()
        if not hasattr(self, "_history"):
            self._history = {}
        self._history[key] = {"action": action,"details": details}

    def backup(self,fullhistory=None,message=None):
        """Backup current state in _previous"""
        previous = self.copy()
        fullhistory = self._fullhistory if fullhistory is None else fullhistory
        if not fullhistory:
            previous._previous = None
        self._previous = previous
        self._events("backup", {"from": self._fullhistory,'message': message})

    def restore(self):
        """Restore the previous signal version if available"""
        if hasattr(self, "_previous") and isinstance(self._previous, signal):
            restored = self._previous.copy()
            for attr in ['x', 'y', 'name', 'type', 'x_label', 'x_unit',
                         'y_label', 'y_unit', 'color', 'linestyle',
                         'linewidth', 'source', 'metadata', '_previous', '_history']:
                setattr(self, attr, getattr(restored, attr))
            self._events("restore", {"from": restored.name})
        else:
            raise AttributeError("No previous signal to restore")

    def enable_fullhistory(self):
        """Enable full history tracking"""
        self._events("enable full history", {"from": "enable_fullhistory"})
        set._fullhistory = True

    def disable_fullhistory(self):
        """Disable full history tracking"""
        set._fullhistory = True
        self._previous = False
        self._events("disable full history", {"from": "disable_fullhistory"})

    @classmethod
    def from_peaks(cls, peaks_obj, x=None, generator_map=None, name="from_peaks", x0=None, n=1000):
        """
        Generate a signal from a set of peaks.

        Parameters
        ----------
        peaks_obj : peaks
            A list-like object containing peak definitions.
        x : array-like, float, or None
            If None: compute x domain from peaks.
            If scalar: interpreted as xmax; linspace from x0 to xmax.
            If array: use as x directly.
        generator_map : dict or None
            Optional map of peak type → generator instance (default is Gaussian).
        name : str
            Name of the signal instance.
        x0 : float or None
            Left bound of the domain (used only if x is None or scalar).
            If None: inferred from peaks.
        n : int
            Number of points in the generated x array.

        Returns
        -------
        signal
            A new signal instance generated from the peaks.

        Example
        -------
        p = peaks()
        p.add(x=[400, 800, 1600], w=30, h=[1.0, 0.6, 0.9], type="gauss")
        s = signal.from_peaks(p, x0=300, n=2048)
        s.plot()
        """
        if x is None:
            xmin = min(p['x'] - 3 * p['w'] for p in peaks_obj)
            xmax = max(p['x'] + 3 * p['w'] for p in peaks_obj)
            if x0 is not None:
                xmin = x0
            x = np.linspace(xmin, xmax, n)
        elif np.isscalar(x):
            xmin = x0 if x0 is not None else min(p['x'] - 3 * p['w'] for p in peaks_obj)
            x = np.linspace(xmin, float(x), n)
        else:
            x = np.asarray(x)

        y = np.zeros_like(x)
        generator_map = generator_map or {}
        for p in peaks_obj:
            g = generator_map.get(p['type'], generator(p['type']))
            y += g(x, p['x'], p['w'], p['h'])
        s = cls(x, y, name=name)
        s.source = "peaks"
        return s

    def sample(self, x_new):
        """Interpolate values from x"""
        return np.interp(x_new, self.x, self.y)

    def plot(self, ax=None, label=None, color=None, linestyle=None, linewidth=None,
             fontsize=12, newfig=False):
        """
        Plot the signal using matplotlib, applying either internal style settings
        or overrides provided at call time.

        Parameters
        ----------
        ax : matplotlib.axes.Axes, optional
            Axis to plot on. If None, uses current axis or new figure if newfig=True.
        label : str, optional
            Legend label. Defaults to self.name.
        color : str or None
            Line color. If None, uses default matplotlib cycling.
        linestyle : str or None
            Line style (e.g., '-', '--'). If None, uses self.linestyle.
        linewidth : float or None
            Line width. If None, uses self.linewidth.
        fontsize : int or str
            Font size for axis labels and legend. Can use values like 'small', 'large'.
        newfig : bool
            If True, creates a new figure before plotting.

        Returns
        -------
        matplotlib.figure.Figure
        matplotlib.axes.Axes
        """
        # Convert font size keywords to numeric
        fontsize_map = {
            "xx-small": 6, "x-small": 8, "small": 10,
            "medium": 12, "large": 14, "x-large": 16, "xx-large": 18
        }
        if isinstance(fontsize, str):
            fontsize = fontsize_map.get(fontsize.lower(), 12)

        if newfig:
            fig = plt.figure()
        else:
            fig = plt.gcf()

        ax = ax or plt.gca()

        # Final style resolution
        label = label if label is not None else self.name
        color = color if color is not None else self.color
        linestyle = linestyle if linestyle is not None else self.linestyle
        linewidth = linewidth if linewidth is not None else self.linewidth

        # Plot
        if self.x is None:
            ax.plot(self.y, label=label, color=color,
                    linestyle=linestyle, linewidth=linewidth)
        else:
            ax.plot(self.x, self.y, label=label, color=color,
                    linestyle=linestyle, linewidth=linewidth)

        # Axes labeling and legend
        if self.x_label and self.x_unit:
            ax.set_xlabel(f"{self.x_label} [{self.x_unit}]", fontsize=fontsize)
        elif self.x_label:
            ax.set_xlabel(self.x_label, fontsize=fontsize)

        if self.y_label and self.y_unit:
            ax.set_ylabel(f"{self.y_label} [{self.y_unit}]", fontsize=fontsize)
        elif self.y_label:
            ax.set_ylabel(self.y_label, fontsize=fontsize)

        if label:
            ax.legend(fontsize=fontsize)

        return fig, ax


    def __repr__(self):
        def fmt_str(s, maxlen=60):
            s = str(s)
            if len(s) <= maxlen:
                return s
            return f"{s[:maxlen//2 - 2]}...{s[-maxlen//2 + 1:]}"
        meta = self.metadata
        span = (f"{self.x[0]:.2f}", f"{self.x[-1]:.2f}") if self.x is not None else ("?", "?")
        size = len(self.x) if self.x is not None else self.n
        field_width = 10
        lines = [
            f"<signal '{self.name}' [{self.type}]>",
            f"{'domain:'.ljust(field_width)} {self.x_label} [{self.x_unit}], span: {span[0]} → {span[1]}, points: {size}",
            f"{'source:'.ljust(field_width)} {self.source}",
            f"{'created:'.ljust(field_width)} {fmt_str(meta.get('date', '?'))}",
            f"{'user:'.ljust(field_width)} {meta.get('user', '?')}@{meta.get('host', '?')}",
            f"{'cwd:'.ljust(field_width)} {fmt_str(meta.get('cwd', '?'))}",
            f"{'version:'.ljust(field_width)} {meta.get('version', '?')}"
        ]
        print("\n".join(lines))
        return str(self)

    def __str__(self):
        if self.y is None:
            return f"<empty {self.type}-signal - source='{self.source}'>"
        else:
            return f"<{self.name}:{self.type}-signal of length={self.n} source='{self.source}'>"

    def align_with(self, other, mode='union', n=1000):
        """
        Align two signals to a common x grid with interpolation and padding.

        Parameters:
            other (signal): the other signal to align with
            mode (str): 'union' (default) or 'intersection'
            n (int): number of points for the new grid

        Returns:
            tuple: (self_interp, other_interp) as new signal instances
        """
        if not isinstance(other, signal):
            raise TypeError("Can only align with another signal.")

        # Determine new x-axis
        if mode == 'union':
            x_min = min(self.x[0], other.x[0])
            x_max = max(self.x[-1], other.x[-1])
        elif mode == 'intersection':
            x_min = max(self.x[0], other.x[0])
            x_max = min(self.x[-1], other.x[-1])
            if x_min >= x_max:
                raise ValueError("No overlapping x-range for 'intersection' mode.")
        else:
            raise ValueError("mode must be 'union' or 'intersection'")

        x_new = np.linspace(x_min, x_max, n)

        # Interpolate and zero outside original domain
        def interp_with_padding(sig):
            y_new = np.interp(x_new, sig.x, sig.y, left=0, right=0)
            return signal(x_new, y_new, name=sig.name, source=sig.source)

        return interp_with_padding(self), interp_with_padding(other)

    def add_noise(self, kind="gaussian", scale=1.0, bias=None):
        """Return a new signal with noise and/or bias added."""
        new_y = self.y.copy()

        # --- Apply bias ---
        if bias is not None:
            if isinstance(bias, (int, float)):
                new_y += bias
            elif isinstance(bias, np.ndarray):
                if bias.shape != self.x.shape:
                    raise ValueError("Bias array must match x domain.")
                new_y += bias
            elif bias == "ramp":
                ramp = np.linspace(0, 1, len(self.x))
                new_y += ramp
            elif isinstance(bias, signal):
                interp_bias = np.interp(self.x, bias.x, bias.y)
                new_y += interp_bias
            else:
                raise TypeError("Invalid bias type.")

        # --- Apply noise ---
        rng = np.random.default_rng()
        if kind == "gaussian":
            new_y += rng.normal(loc=0.0, scale=scale, size=self.y.shape)
        elif kind == "poisson":
            if np.any(new_y < 0):
                raise ValueError("Poisson noise requires non-negative values.")
            new_y = rng.poisson(lam=new_y * scale) / scale
        else:
            raise ValueError("Unknown noise type.")

        return signal(self.x.copy(), new_y, name=self.name + "+noise", source=self.source + "+noise")

    def copy(self):
        """Deep copy of the signal, excluding full history control flag"""
        new = signal(
            x=self.x.copy() if self.x is not None else None,
            y=self.y.copy() if self.y is not None else None,
            name=self.name,
            type=self.type,
            x_label=self.x_label,
            x_unit=self.x_unit,
            y_label=self.y_label,
            y_unit=self.y_unit,
            metadata=self.metadata.copy(),
            source=self.source,
            color=self.color,
            linewidth=self.linewidth,
            linestyle=self.linestyle,
            fullhistory=self._fullhistory
        )
        new._history = self._history.copy()
        if self._fullhistory and self._previous is not None:
            new._previous = self._previous.copy()
        else:
            new._previous = None
        return new

    def _binary_op(self, other, op):
        """Binary operation on signals"""
        if isinstance(other, (int, float)):
            return signal(self.x.copy(), op(self.y, other), name=self.name, source=self.source)

        if isinstance(other, signal):
            s1, s2 = self.align_with(other)
            return signal(s1.x, op(s1.y, s2.y), name=f"({s1.name}){op.__name__}({s2.name})")

        raise TypeError(f"Unsupported operand type(s) for {op.__name__}: 'signal' and '{type(other).__name__}'")

    def __add__(self, other): return self._binary_op(other, operator.add)
    def __sub__(self, other): return self._binary_op(other, operator.sub)
    def __mul__(self, other): return self._binary_op(other, operator.mul)
    def __truediv__(self, other): return self._binary_op(other, operator.truediv)
    def __iadd__(self, other): return self._binary_op(other, operator.add)
    def __isub__(self, other): return self._binary_op(other, operator.sub)
    def __imul__(self, other): return self._binary_op(other, operator.mul)
    def __itruediv__(self, other): return self._binary_op(other, operator.truediv)

    def _to_serializable(self):
        """Convert the signal into a dictionary suitable for JSON export."""
        return {
            "x": self.x.tolist() if self.x is not None else None,
            "y": self.y.tolist() if self.y is not None else None,
            "name": self.name,
            "type": self.type,
            "x_label": self.x_label,
            "x_unit": self.x_unit,
            "y_label": self.y_label,
            "y_unit": self.y_unit,
            "source": self.source,
            "metadata": self.metadata,
            "color": self.color,
            "linestyle": self.linestyle,
            "linewidth": self.linewidth,
            "_history": self._history,
            "_previous": self._previous.to_serializable() if self._fullhistory and self._previous else None
        }

    @staticmethod
    def _from_serializable(data,message=None):
        """Convert a serialized dict to a signal"""
        s = signal(
            x=np.array(data["x"]) if data["x"] is not None else None,
            y=np.array(data["y"]) if data["y"] is not None else None,
            name=data.get("name", "signal"),
            type=data.get("type", "generic"),
            x_label=data.get("x_label", "index"),
            x_unit=data.get("x_unit", "-"),
            y_label=data.get("y_label", "intensity"),
            y_unit=data.get("y_unit", "a.u."),
            metadata=data.get("metadata", {}),
            source=data.get("source", "array"),
            color=data.get("color", None),
            linestyle=data.get("linestyle", "-"),
            linewidth=data.get("linewidth", 2),
            fullhistory=True,  # Set default on load
            message = message
        )
        s._history = data.get("_history", {})
        prev = data.get("_previous", None)
        s._previous = signal._from_serializable(prev) if prev else None
        return s

    def save(self, filepath=None, zip=True, export_csv=False):
        """
        Save signal to JSON (optionally compressed) and optionally CSV.

        Parameters
        ----------
        filepath : str or Path or None
            If None, builds path from metadata['cwd'] and self.name + '.json[.gz]'.
            If a directory, appends name + '.json[.gz]'.
            If a file, uses as is.
        zip : bool
            Whether to compress the JSON file using gzip. Default: True.
        export_csv : bool
            If True, also save a .csv file (x,y) alongside the JSON.
        """
        # Resolve filepath
        if filepath is None:
            filepath = os.path.join(self.metadata["cwd"], f"{self.name}.json.gz" if zip else f"{self.name}.json")
        elif os.path.isdir(filepath):
            filepath = os.path.join(filepath, f"{self.name}.json.gz" if zip else f"{self.name}.json")
        filepath = Path(filepath)
        data = self._to_serializable()
        if zip or str(filepath).endswith(".gz"):
            with gzip.open(filepath, "wt", encoding="utf-8") as f:
                json.dump(data, f, indent=2)
        else:
            with open(filepath, "w", encoding="utf-8") as f:
                json.dump(data, f, indent=2)

        # Optionally export to CSV
        if export_csv and self.x is not None and self.y is not None:
            csv_path = filepath.with_suffix(".csv")
            with open(csv_path, "w", encoding="utf-8") as f:
                f.write("x,y\n")
                for xi, yi in zip(self.x, self.y):
                    f.write(f"{xi},{yi}\n")

    @staticmethod
    def load(filepath):
        """
        Load a signal from a JSON or gzipped JSON file, including recursive _previous.

        Parameters
        ----------
        filepath : str or Path
            Path to the JSON or .gz file

        Returns
        -------
        signal
            A fully reconstructed signal object
        """
        filepath = Path(filepath)
        if filepath.suffix == ".gz":
            with gzip.open(filepath, "rt", encoding="utf-8") as f:
                data = json.load(f)
        else:
            with open(filepath, "r", encoding="utf-8") as f:
                data = json.load(f)
        return signal._from_serializable(data,message=f"loaded from {filepath}")

    def apply_poisson_baseline_filter(self, window_ratio=0.02, gain=1.0, proba=0.9):
        """
        Apply a baseline filter assuming Poisson-dominated statistics with adjustable gain
        and a rejection threshold based on the Bienaymé-Tchebychev inequality.

        The signal is filtered by removing values likely caused by statistical noise
        (false peaks) using a per-point threshold defined from local statistics:

        - Local mean: $$ \mu_t = \frac{1}{w} \sum_{i \in W(t)} y_i $$
        - Local std dev: $$ \sigma_t = \sqrt{\mu_t \cdot \text{gain}} $$
        - Coefficient of variation: $$ \text{cv}_t = \frac{\sigma_t}{\mu_t} $$
        - Estimated local intensity (lambda): $$ \lambda_t = \frac{1}{\text{cv}_t^2} $$
        - Bienaymé-Tchebychev threshold: $$ \text{threshold}_t = \frac{1}{\sqrt{1 - p}} \cdot \sqrt{10 \lambda_t \cdot \Delta t} $$

        Parameters
        ----------
        window_ratio : float, default=0.02
            Ratio of signal length used as window size (must yield odd integer ≥ 11).
        gain : float, default=1.0
            Linear amplification factor applied to simulate signal counts.
        proba : float, default=0.9
            Minimum probability to consider a signal point significant. Must be in (0, 1).

        Returns
        -------
        signal
            The current signal instance (self), with updated `y`.

        Raises
        ------
        ValueError
            If the window size is too small for reliable statistics.
        """
        import numpy as np

        y = self.y
        n = len(y)

        # --- Window size ---
        w = int(window_ratio * n)
        if w < 11:
            raise ValueError(f"Window too small ({w} < 11); increase signal length or window_ratio.")
        if w % 2 == 0:
            w += 1

        # --- Sliding window views ---
        padded = np.pad(y, w//2, mode='reflect')
        windows = sliding_window_view(padded, w)  # shape: (n, w)
        # --- Robust local statistics ---
        local_mean = np.mean(windows, axis=1)
        local_std = np.std(windows, axis=1)
        # Prevent division by zero
        eps = 1e-12
        cv = np.where(local_mean > eps, local_std / local_mean, np.inf)
        # Estimate lambda from CV (avoid inf/NaN)
        lam = np.where(cv > 0, 1 / (cv ** 2), 0)
        # Bienaymé-Tchebychev threshold with gain
        delta_t = self.sampling_dt if hasattr(self, "sampling_dt") else 1.0
        k = 1 / np.sqrt(1 - proba)
        threshold = k * np.sqrt(lam * delta_t * gain)

        # --- Backup and History ---
        self.backup()
        event = f"apply_poisson_baseline_filter(window={w}, proba={proba:.2f}, gain={gain})"
        self._events("filter", {"from": event})

        # --- Apply filter ---
        self.y = np.where(self.y > threshold, self.y, 0.0)
        return self


    def _toDNA(self,encode=True,scales=[1,3,4,8,16,32]):
        """
        Return a DNA encoded signal

        Parameters
            encode : bool (default=True)
            scales : list (default=[1,3,4,8,16,32])
        """
        return DNAsignal(self,encode=encode)

# ------------------------
# Signal collection class
# ------------------------
class signal_collection(list):
    """
    A container class for multiple `signal` instances that ensures alignment on a shared x-grid.

    The collection is used to manage, compare, combine, or visualize multiple signals (e.g., from
    replicates, experiments, synthetic scenarios). Signals are interpolated and padded on insertion
    so all have the same shape and domain. Arithmetic, matrix extraction, and overlay plots are supported.

    Parameters:
    -----------
    *signals : signal
        One or more signal instances to include (they are copied and aligned).
    n : int
        Number of sampling points in the aligned x grid (default: 1000).
    mode : str
        Alignment mode: 'union' or 'intersection' of x-ranges.

    Core Attributes:
    ----------------
    mode : str
        Alignment strategy used ("union" or "intersection").
    n : int
        Number of x-points used in alignment (default=1024).

    Key Methods:
    ------------
    - append(signal)              → add and align a new signal
    - to_matrix()                 → convert signals to a 2D array (n_signals x n_points)
    - mean(coeffs=None)           → weighted or unweighted mean
    - sum(coeffs=None)            → weighted or unweighted sum
    - plot(...)                   → overlay signals with optional mean/sum
    - copy                        → all signals stored are deep copies
    - generate_synthetic          → signal collection composed of random peaks.
    - __getitem__(...)            → slice, list, or name-based access to signals
    - __repr__ / __str__          → report contents with span and names
    - _toDNA(signal_collection)   → list of DNAsignals

    Access Patterns:
    ----------------
    - sc[0:3]           → subcollection by slice
    - sc[[0, 2]]        → subcollection by list of indices
    - sc["name"]        → return a copy of signal with that name
    - sc["A", "B"]      → return a subcollection with those names

    Supports arithmetic operations for aligned signal mixtures.

    Arithmetic operations on aligned signal collections
    ----------------------------------------------------
    - Scalar multiplication: `a * sc` scales each signal by a constant `a`.
    - Collection addition: `sc1 + sc2` adds two collections element-wise.
    - Linear combinations: `a * A + b * B + c * C` constructs mixtures of compatible collections.
    - Compatible with `sum([a*A, b*B, ...])` for aggregating multiple weighted collections.

    Requirements
        - All signal_collections must share the same number of signals.
        - Signals are aligned on a common x-grid (same `n`, `mode`, and domain).
        - Element-wise operations preserve signal names and metadata when possible.


    Examples:
    ---------
    >>> sc = signal_collection(s1, s2, s3)
    >>> sc.plot(show_mean=True)

    >>> sc[0:2]         # sub-collection (copy)
    >>> sc["peak1"]     # get copy of signal named 'peak1'
    >>> mat = sc.to_matrix()

    >>> sc.mean().plot()
    >>> sc.sum(coeffs=[0.4, 0.6]).plot()
    """

    def __init__(self, *signals, n=1024, mode='union', name=None, force=True):
        """
        Initialize collection with aligned signals of the same type.

        Parameters
        ----------
        *signals : signal
            Signal instances to include.
        n : int (default=1000)
            Number of points on the common x grid.
        mode : str (default="union")
            'union' or 'intersection' for domain alignment.
        force : bool (default=True)
            If False, require all signals to have the same type.
        """
        nsignals = len(signals)
        self.name = f"Collection of {nsignals} signals" if name is None else name
        self.mode = mode
        self.n = n
        if not force and nsignals > 1:
            types = {s.type for s in signals}
            if len(types) > 1:
                raise ValueError(f"Incompatible signal types: {types}. Use force=True to override.")
        s_copy = [s.copy() for s in signals]
        for s in s_copy:
            if s.x is None:
                s.x = np.linspace(0,s.n-1,s.n,endpoint=True,dtype=np.uint32)
        aligned = self._align_all(s_copy)
        super().__init__(aligned)
        self._plotted_once = False # flag to track plotting

    def append(self, new_signal):
        """Append and align the new signal to the existing collection."""
        aligned = self._align_all(list.__add__(self, [new_signal.copy()]))
        self.clear()
        self.extend(aligned)

    def __getitem__(self, key):
        """
            sc["my_signal"]                  # returns a copy of the signal named "my_signal"
            sc["A", "B", "C"]                # returns a sub-collection with those names
            sc[0:2] or sc[[0, 2]]            # still works for index-based access
        """
        if isinstance(key, slice):
            return signal_collection(*(s.copy() for s in super().__getitem__(key)))
        elif isinstance(key, list):
            return signal_collection(*(self[k].copy() for k in key))
        elif isinstance(key, tuple):
            return signal_collection(*(s.copy() for s in self if s.name in key))
        elif isinstance(key, str):
            for s in self:
                if s.name == key:
                    return s.copy()
            raise KeyError(f"No signal named '{key}' in collection.")
        else:
            return super().__getitem__(key)


    def __setitem__(self, key, value):
        if not isinstance(value, signal):
            raise TypeError("Only signal instances can be assigned.")
        aligned = self._align_all(self[:key] + [value.copy()] + self[key+1:])
        super().__setitem__(key, aligned[key])

    def __delitem__(self, key):
        if isinstance(key, list):
            # delete in reverse to preserve indices
            for k in sorted(key, reverse=True):
                super().__delitem__(k)
        else:
            super().__delitem__(key)

    def _align_all(self, signals):
        if not signals:
            return []
        # Determine common grid
        x_min = min(s.x[0] for s in signals)
        x_max = max(s.x[-1] for s in signals)
        x_common = np.linspace(x_min, x_max, self.n)
        # Interpolate all to shared grid
        aligned = [signal(x_common, np.interp(x_common, s.x, s.y, left=0, right=0),
                          name=s.name, source=s.source) for s in signals]
        return aligned

    def __str__(self):
        return f"<signal_collection with {len(self)} aligned signals>"

    def __repr__(self):
        lines = [f"[{i}] '{s.name}' span=({s.x[0]:.2f}, {s.x[-1]:.2f})" for i, s in enumerate(self)]
        print("<signal_collection>\n" + "\n".join(lines))
        return str(self)

    def to_matrix(self):
        """Return a 2D array (n_signals x n_points) of aligned signal values."""
        return np.vstack([s.y for s in self])

    def sum(self, indices_or_names=None, coeffs=None):
        """
        Sum selected signals, optionally weighted by coeffs.

        Parameters
        ----------
        indices_or_names : list[int or str], optional
            If provided, selects a subset by index or name.
        coeffs : list[float], optional
            Weights matching the number of selected signals.

        Returns
        -------
        signal
            Summed signal.
        """
        if indices_or_names is None:
            selected = self
        else:
            selected = []
            for k in indices_or_names:
                if isinstance(k, int):
                    selected.append(self[k])
                elif isinstance(k, str):
                    found = next((s for s in self if s.name == k), None)
                    if found is None:
                        raise KeyError(f"Signal name '{k}' not found.")
                    selected.append(found)
                else:
                    raise TypeError("indices_or_names must contain ints or strs.")

        matrix = np.vstack([s.y for s in selected])
        if coeffs is not None:
            coeffs = np.asarray(coeffs)
            if coeffs.shape[0] != len(selected):
                raise ValueError("Length of coeffs must match number of selected signals.")
            result = np.dot(coeffs, matrix)
        else:
            result = np.sum(matrix, axis=0)

        return signal(selected[0].x.copy(), result, name="sum", source="signal_collection")

    def mean(self, indices_or_names=None, coeffs=None):
        """
        Mean of selected signals, optionally weighted.

        Parameters
        ----------
        indices_or_names : list[int or str], optional
            Signal names or indices to include.
        coeffs : list[float], optional
            Weights for selected signals.

        Returns
        -------
        signal
            Averaged signal.
        """
        if coeffs is not None:
            total_weight = np.sum(coeffs)
        else:
            total_weight = len(indices_or_names) if indices_or_names else len(self)

        s = self.sum(indices_or_names=indices_or_names, coeffs=coeffs)
        return signal(s.x.copy(), s.y / total_weight, name="mean", source="signal_collection")

    def plot(self, indices=None, labels=True, title=None, newfig=None, ax=None,
             show_mean=False, show_sum=False, coeffs=None, fontsize=12, colormap=None):
        """
        Plot selected signals with style attributes and optional overlays.

        Parameters
        ----------
        indices : list[int] or list[str], optional
            Signals to plot by index or name.
        labels : bool
            Whether to show signal labels.
        title : str
            Plot title.
        newfig : bool or None
            If True, always open a new figure.
            If False, use current axes.
            If None, open new figure only the first time this collection is plotted.
        ax : matplotlib axis, optional
            Axis to draw on.
        show_mean : bool
            Overlay mean curve.
        show_sum : bool
            Overlay sum curve.
        coeffs : list[float], optional
            Optional weights for mean/sum.
        fontsize : int or str
            Font size for labels and legend.
        colormap : list[str], optional
            List of colors to cycle through when signal.color is None.

        Returns
        -------
        matplotlib.figure.Figure
        """
        if not self:
            return
        # Select signals to plot
        if indices is None:
            selected = self
        else:
            if isinstance(indices[0], str):
                selected = [s for s in self if s.name in indices]
            else:
                selected = [self[i] for i in indices]
        # Circular color iterator
        if colormap is None:
            colormap = plt.rcParams['axes.prop_cycle'].by_key().get('color', [])
        color_cycle = (colormap * ((len(selected) // len(colormap)) + 1)) if colormap else [None] * len(selected)
        # Determine whether to create new figure
        if newfig is None:
            newfig = not self._plotted_once
        if newfig:
            fig = plt.figure()
        else:
            fig = plt.gcf()
        self._plotted_once = True
        ax = ax or plt.gca()
        for i, s in enumerate(selected):
            c = s.color if s.color is not None else color_cycle[i]
            s.plot(ax=ax, label=s.name if labels else None,
                   color=c, linestyle=s.linestyle, linewidth=s.linewidth, fontsize=fontsize)
        # Overlay mean and sum
        if show_mean:
            self.mean(coeffs).plot(ax=ax, label="mean", color="black", linewidth=2, fontsize=fontsize)
        if show_sum:
            self.sum(coeffs).plot(ax=ax, label="sum", color="red", linestyle="--", linewidth=2, fontsize=fontsize)

        title = self.name if title is None else title
        title = "Signal Collection" if title is None else title
        ax.set_title(title, fontsize=fontsize)
        ax.tick_params(labelsize=fontsize)
        if labels:
            ax.legend(fontsize=fontsize)
        plt.show()
        return fig

    @classmethod
    def generate_synthetic(cls,
                           n_signals=5,
                           n_peaks=5,
                           kind_distribution="uniform",  # or "random"
                           kinds=("gauss", "lorentz", "triangle"),
                           x_range=(0, 1000),
                           n_points=1024,
                           avoid_overlap=True,
                           width_range=(20, 60),
                           height_range=(0.5, 1.0),
                           normalize=True,
                           noise=None,    # e.g. {"kind": "gaussian", "scale": 0.05}
                           bias=None,     # can be scalar, "ramp", or signal instance
                           name_prefix="synthetic",
                           seed=None):
        """
        Generate a synthetic signal collection composed of random peaks.

        Parameters
        ----------
        n_signals : int
            Number of synthetic signals to generate.
        n_peaks : int or tuple[int,int]
            Number of peaks per signal or its range.
        kind_distribution : str
            'uniform' → use all peak kinds equally; 'random' → random draw from kinds.
        kinds : tuple[str]
            Generator types to choose from ('gauss', 'lorentz', 'triangle').
        x_range : tuple[float, float]
            Start and end of x-domain.
        n_points : int
            Number of sampling points for each signal (default: 1024).
        avoid_overlap : bool
            Prevent peaks from overlapping by checking spacing vs. width.
        width_range : tuple[float, float]
            Range of widths for the peaks.
        height_range : tuple[float, float]
            Range of peak heights.
        normalize : bool
            Normalize each signal so the highest peak has intensity 1.
        noise : dict or None
            Optional noise model, e.g. {"kind": "gaussian", "scale": 0.01}.
        bias : float, str, np.ndarray, or signal
            Optional signal bias: can be a constant, 'ramp', or signal.
        name_prefix : str
            Base name for each generated signal.
        seed : int or None
            Random seed for reproducibility.

        Returns
        -------
        signal_collection
            A collection of generated signals.

        Examples
        ---------
        # 1. Default random peaks, Gaussian + ramp bias
        sc = signal_collection.generate_synthetic(
            n_signals=5,
            n_peaks=6,
            kinds=("gauss", "lorentz", "triangle"),
            noise={"kind": "gaussian", "scale": 0.02},
            bias="ramp",
            name_prefix="test"
        )
        sc.plot(show_mean=True, fontsize="large")

        # 2. High-res signal, fixed width and height
        sc2 = signal_collection.generate_synthetic(
            n_signals=3,
            n_peaks=8,
            kinds=("gauss",),
            width_range=(30, 30),
            height_range=(1.0, 1.0),
            x_range=(0, 500),
            n_points=2048,
            normalize=False,
            seed=123
        )
        sc2.plot(fontsize=14)

        # 3. Poisson noise, no overlap, save output
        sc3 = signal_collection.generate_synthetic(
            n_signals=2,
            n_peaks=5,
            noise={"kind": "poisson", "scale": 2.0},
            name_prefix="poisson_example"
        )
        for s in sc3:
            s.save(export_csv=True)
        """
        rng = np.random.default_rng(seed)
        xmin, xmax = x_range
        x = np.linspace(xmin, xmax, n_points)
        all_peak_centers = []
        all_peak_widths = []
        signals = []
        #cover n_peaks into range
        if isinstance(n_peaks,int):
            n_peaks = (n_peaks,n_peaks)

        for i in range(n_signals):
            p = peaks()
            attempts = 0
            current_npeaks = rng.integers(*n_peaks,endpoint=True)
            while len(p) < current_npeaks and attempts < 100 * current_npeaks:
                w = rng.uniform(*width_range)
                h = rng.uniform(*height_range)
                x0 = rng.uniform(xmin + 3 * w, xmax - 3 * w)

                # Check for global overlap
                if avoid_overlap:
                    if any(abs(x0 - xj) < 3 * max(w, wj) for xj, wj in zip(all_peak_centers, all_peak_widths)):
                        attempts += 1
                        continue

                if kind_distribution == "uniform":
                    kind = kinds[len(p) % len(kinds)]
                else:
                    kind = rng.choice(kinds)

                p.add(x=x0, w=w, h=h, type=kind)
                all_peak_centers.append(x0)
                all_peak_widths.append(w)
                attempts += 1

            # warning if the number of peaks is lower than expected
            if len(p) < current_npeaks:
                warnings.warn(f"Signal {i+1}: only placed {len(p)} of {current_npeaks} peaks due to overlap constraints.")

            # Generate signal
            generator_map = {k: generator(k) for k in kinds}
            s = signal.from_peaks(p, x=x, generator_map=generator_map, name=f"{name_prefix}_{i+1}")
            s.y_label = "intensity"
            s.y_unit = "a.u."

            # Normalize
            if normalize:
                ymax = np.max(np.abs(s.y))
                if ymax > 0:
                    s.y /= ymax

            # Add noise or bias
            if noise:
                s = s.add_noise(**noise, bias=bias)
            elif bias is not None:
                s = s.add_noise(kind="gaussian", scale=0.0, bias=bias)  # bias only

            s._source_peaks = p
            signals.append(s)

        return cls(*signals, n=n_points, mode='union', force=True), [s._source_peaks for s in signals]

    @classmethod
    def generate_mixtures(cls,
                          n_mixtures=10,             # int
                          max_peaks = 16,            # int,
                          peaks_per_mixture = (3,8), # tuple[int, int],
                          amplitude_range = (0.5,2), # tuple[float, float],
                          flatten = 'mean',          # literal['sum', 'mean']
                          # parameters
                          n_signals=1,
                          n_peaks=1,
                          kinds=("gauss",),
                          width_range=(0.5, 3),
                          height_range=(1.0, 5.0),
                          x_range=(0, 500),
                          n_points=1024,
                          normalize=False,
                          seed=None,
                          **kwargs):
        """
        Generate synthetic mixtures of signals by combining a subset of base peaks.

        Parameters
        ----------
        n_mixtures : int
            Number of synthetic mixtures to generate.

        max_peaks : int
            Maximum number of base signals (from which peaks are taken).

        peaks_per_mixture : tuple of (int, int)
            Range (min, max) for the number of peaks to combine in each mixture.
            Cannot exceed `max_peaks`.

        amplitude_range : tuple of (float, float)
            Random scaling range applied to peak amplitudes in each mixture.

        flatten : {'sum', 'mean'}, default='mean'
            How to combine the signals for each mixture.

        **kwargs : dict
            All other keyword arguments passed to `generate_synthetic`.

        Returns
        -------
        result_collection : signal_collection
            A collection of synthetic mixed signals.

        all_peaks : list of dict
            All individual peaks originally generated.

        used_peak_ids : list of list of str
            For each mixture, the list of peak names used.

        Examples
        --------
        S, pS = signal_collection.generate_mixtures(
        ...     n_mixtures=30,
        ...     max_peaks=12,
        ...     peaks_per_mixture=(4, 8),
        ...     amplitude_range=(0.2, 1.5),
        ...     n_signals=12,
        ...     kinds=("gauss",),
        ...     width_range=(0.5, 3),
        ...     height_range=(1.0, 5.0),
        ...     x_range=(0, 500),
        ...     n_points=2048,
        ...     normalize=False,
        ...     seed=123
        ... )
        >>> S.plot()
        """

        # Step 1: Set defaults and validate
        if not (0 < peaks_per_mixture[0] <= peaks_per_mixture[1] <= max_peaks):
            raise ValueError("Invalid peaks_per_mixture range.")

        if seed is not None:
            random.seed(seed)
            np.random.seed(seed)

        # Step 2: Generate base signals and flatten peaks
        S, peak_list = cls.generate_synthetic(
            n_signals=max_peaks,
            n_peaks=n_peaks,
            kinds=("gauss",),
            width_range=width_range,
            height_range=height_range,
            x_range=x_range,
            n_points=n_points,
            normalize=normalize,
            seed=seed,
            **kwargs
            )
        peak_objs = peaks(peak_list)  # Flatten into a peak collection

        all_peaks = list(peak_objs)
        peak_id_map = {p['name']: p for p in all_peaks}

        mixtures = []
        used_peak_ids = []

        for i in range(n_mixtures):
            n_peaks_in_mix = random.randint(*peaks_per_mixture)
            selected_peaks = random.sample(all_peaks, n_peaks_in_mix)
            mixture_signal = np.zeros_like(S[0].y)

            current_ids = []
            for p in selected_peaks:
                scale = random.uniform(*amplitude_range)
                yvals = generator(p['type'])(S[0].x, p['x'], p['w'], p['h'] * scale)
                mixture_signal += yvals
                current_ids.append(p['name'])

            if flatten == 'mean':
                mixture_signal /= len(selected_peaks)

            sig = signal(x=S[0].x, y=mixture_signal, name=f"mixture {i} of {n_mixtures}")
            mixtures.append(sig)
            used_peak_ids.append(current_ids)

        result_collection = cls(*mixtures,n=n_points, mode='union', force=True)

        return result_collection, all_peaks, used_peak_ids


    def _toDNA(self,encode=True,scales=[1,3,4,8,16,32]):
        """
        Return a collection of DNA encoded signals  (class DNAsignal_collection)

        Parameters
            encode : bool (default=True)
            scales : list (default=[1,3,4,8,16,32])
        """
        return DNAsignal_collection(*[DNAsignal(s,scales=scales,encode=encode) for s in self])

    def __mul__(self, scalar):
        """override *"""
        if not isinstance(scalar, (int, float)):
            raise TypeError("Only scalar multiplication is supported.")
        new_signals = [signal(s.x, s.y * scalar, name=s.name, source=s.source) for s in self]
        return signal_collection(*new_signals, n=self.n, mode=self.mode, name=f"{scalar}*{self.name}")

    __rmul__ = __mul__  # allow scalar multiplication from left

    def __add__(self, other):
        """overrride +"""
        if not isinstance(other, signal_collection):
            raise TypeError("Can only add another signal_collection.")
        if len(self) != len(other):
            raise ValueError("Signal collections must have the same number of signals.")
        new_signals = []
        for s1, s2 in zip(self, other):
            if not np.allclose(s1.x, s2.x):
                raise ValueError(f"Signals '{s1.name}' and '{s2.name}' have incompatible x-grids.")
            new_y = s1.y + s2.y
            new_name = f"{s1.name}+{s2.name}"
            new_signals.append(signal(s1.x.copy(), new_y, name=new_name, source="add"))
        return signal_collection(*new_signals, n=self.n, mode=self.mode, name=f"({self.name})+({other.name})")

    def __radd__(self, other):
        """override +"""
        if other == 0:  # allow sum([...]) to work
            return self
        return self.__add__(other)

    def __iadd__(self, other):
        """override +="""
        return self.__add__(other)

# %% Peaks class
# ------------------------
class peaks:
    """
    A class for managing a collection of peak definitions used in synthetic signal generation.

    Each peak is represented as a dictionary with the following fields:
    - 'name' (str): unique identifier (autogenerated if not provided)
    - 'x' (float): center position (e.g., time, wavenumber, index)
    - 'w' (float): width (related to FWHM)
    - 'h' (float): peak height
    - 'type' (str): generator type (e.g., 'gauss', 'lorentz', 'triangle')

    Supports:
    - Flexible addition and broadcasting of peak parameters
    - Named or indexed access to individual or multiple peaks
    - Overloaded operators for peak translation and scaling
    - Utility methods: update, sort, rename, remove_duplicates, copy
    - Conversion to signal object via `.to_signal()`
    - Informative __str__ and __repr__ output

    This class is used to build reproducible and structured test cases for symbolic encoding (e.g., sig2dna).
    """

    def __init__(self, data=None):
        """
        Initialize a peak collection.

        Parameters
        ----------
        data : list of dict, list of peaks, or None
            - If None: create an empty peaks object.
            - If list of dict: must contain at least 'x', 'w', 'h'; 'name' and 'type' are optional.
            - If list of peaks: flattens into a single collection.

            If provided, must be a list of dictionaries with the keys:
            - 'name' : str (optional; auto-generated if missing or duplicated)
            - 'x'    : float (center position)
            - 'w'    : float (width)
            - 'h'    : float (height)
            - 'type' : str (e.g., 'gauss')

        If data is None, initializes an empty peaks object.
        """
        self._peaks = []
        self._names = set()
        if data:
            flat_list = []
            if all(isinstance(p, peaks) for p in data):
                # Merge from list of peaks instances
                for p_obj in data:
                    flat_list.extend(p_obj._peaks)
            else:
                flat_list = data
            for p in flat_list:
                name = p.get("name", f"P{len(self._peaks)}")
                while name in self._names:
                    base = name or f"P{len(self._peaks)}"
                    suffix = 0
                    while True:
                        unique_name = f"{base}_{suffix}" if suffix > 0 else base
                        if unique_name not in self._names:
                            break
                        suffix += 1
                    name = unique_name
                peak = {
                    "name": name,
                    "x": float(p["x"]),
                    "w": float(p["w"]),
                    "h": float(p["h"]),
                    "type": p.get("type", "gauss")
                }
                self._peaks.append(peak)
                self._names.add(name)
                self.sort()

    def update(self, data):
        """
        Update or insert peaks from a list of dictionaries.

        Parameters
        ----------
        data : list of dict
            Each dict must include at least 'x', 'w', 'h'.
            If 'name' matches an existing peak, it will be updated.
            If 'name' is new or missing, the peak is appended.
        """
        for p in data:
            name = p.get("name", None)
            peak = {
                "name": name,
                "x": float(p["x"]),
                "w": float(p["w"]),
                "h": float(p["h"]),
                "type": p.get("type", "gauss")
            }

            if name and name in self._names:
                for i, existing in enumerate(self._peaks):
                    if existing["name"] == name:
                        self._peaks[i].update(peak)
                        break
            else:
                peak["name"] = name or f"P{len(self._peaks)}"
                while peak["name"] in self._names:
                    peak["name"] += "_"
                self._peaks.append(peak)
                self._names.add(peak["name"])
        return self

    def add(self, x, w=1.0, h=1.0, name=None, type="gauss"):
        """
        Add one or multiple peaks to the collection.

        Parameters
        ----------
        x : float or array-like
            Center positions of the peaks.
        w : float or array-like
            Width(s) of the peaks (broadcastable).
        h : float or array-like
            Height(s) of the peaks (broadcastable).
        name : str or list of str or None
            Peak name(s); auto-generated if None or duplicate.
        type : str
            Generator type, e.g., 'gauss', 'lorentz', etc.
        """
        x = np.atleast_1d(x)
        w = np.full_like(x, w) if np.isscalar(w) else np.asarray(w)
        h = np.full_like(x, h) if np.isscalar(h) else np.asarray(h)
        names = name if isinstance(name, (list, tuple)) else [name] * len(x)

        for i, (xi, wi, hi, ni) in enumerate(zip(x, w, h, names)):
            base = ni or f"P{len(self._peaks)}"
            suffix = 0
            while True:
                candidate = f"{base}_{suffix}" if suffix > 0 else base
                if candidate not in self._names:
                    break
                suffix += 1
            ni = candidate

            self._peaks.append({
                "name": ni,
                "x": float(xi),
                "w": float(wi),
                "h": float(hi),
                "type": type
            })
            self._names.add(ni)

    def __getitem__(self, key):
        if isinstance(key, str):
            next((p for p in self._peaks if p['name'] == key), None)
        elif isinstance(key, int):
            return self._peaks[key]
        elif isinstance(key, (list, tuple)):
            return [self[k] for k in key]
        elif isinstance(key, slice):
            return self._peaks[key]
        else:
            raise KeyError("Unsupported index type")

    def __delitem__(self, key):
        if isinstance(key, str):
            self._peaks = [p for p in self._peaks if p['name'] != key]
            self._names.discard(key)
        elif isinstance(key, int):
            self._names.discard(self._peaks[key]['name'])
            del self._peaks[key]
        elif isinstance(key, (list, tuple)):
            for k in key:
                self.__delitem__(k)
        elif isinstance(key, slice):
            for p in self._peaks[key]:
                self._names.discard(p['name'])
            del self._peaks[key]

    def __add__(self, shift):
        p = peaks()
        for peak in self._peaks:
            new_peak = peak.copy()
            new_peak['x'] += shift if np.isscalar(shift) else shift.pop(0)
            p._peaks.append(new_peak)
            p._names.add(new_peak['name'])
        return p

    def __mul__(self, factor):
        p = peaks()
        for peak in self._peaks:
            new_peak = peak.copy()
            if isinstance(factor, tuple):
                new_peak['w'] *= factor[0]
                new_peak['h'] *= factor[1]
            else:
                f = factor if np.isscalar(factor) else factor.pop(0)
                new_peak['w'] *= f
                new_peak['h'] *= f
            p._peaks.append(new_peak)
            p._names.add(new_peak['name'])
        return p

    def __truediv__(self, factor):
        inv = (1/factor[0], 1/factor[1]) if isinstance(factor, tuple) else 1/factor
        return self * inv

    def rename(self, prefix="P"):
        self._names.clear()
        for i, peak in enumerate(self._peaks):
            peak['name'] = f"{prefix}{i}"
            self._names.add(peak['name'])

    def remove_duplicates(self):
        seen = set()
        unique_peaks = []
        for peak in self._peaks:
            key = (peak['x'], peak['w'], peak['h'], peak['type'])
            if key not in seen:
                seen.add(key)
                unique_peaks.append(peak)
        self._peaks = unique_peaks
        self._names = {p['name'] for p in self._peaks}

    def __len__(self):
        return len(self._peaks)

    def __repr__(self):
        if not self._peaks:
            return str(self)

        # Compute max name width
        name_width = max(len(p['name']) for p in self._peaks)
        type_width = max(len(p['type']) for p in self._peaks)

        # Format header (optional)
        header = f"{'name':<{name_width}}  {'x':>8}  {'w':>8}  {'h':>8}  {'type':<{type_width}}"
        lines = [header, "-" * len(header)]

        # Format each peak
        for p in self._peaks:
            lines.append(
                f"{p['name']:<{name_width}}  "
                f"{p['x']:>8.2f}  {p['w']:>8.2f}  {p['h']:>8.2f}  {p['type']:<{type_width}}"
            )

        return "\n".join(lines)

    def __str__(self):
        if not self._peaks:
            return "<peaks instance with 0 peaks>"
        xmin = min(p['x'] for p in self._peaks)
        xmax = max(p['x'] for p in self._peaks)
        return f"<peaks instance with {len(self._peaks)} peaks> spanned from x={xmin:.2f} to x={xmax:.2f}"

    def names(self):
        """Return the list of names"""
        return [p['name'] for p in self._peaks]

    def as_dict(self):
        """Return the list of peaks as dict"""
        return self._peaks.copy()

    def to_signal(self, index=None, name=None, generator_map=None, x=None, x0=0.0, n=1000):
        """Generate a signal from a peaks object. Optionally restrict to a subset."""
        selected = self._peaks
        if index is not None:
            if isinstance(index, (str, int)):
                index = [index]
            selected = [self[k] for k in index]
        return signal.from_peaks(selected, x=x, name=name or "from_peaks", generator_map=generator_map, x0=x0, n=n)

    def copy(self):
        """Return a deep-copy of the peaks"""
        new = peaks()
        new._peaks = deepcopy(self._peaks)
        return new

    def sort(self, order="asc"):
        """
        Sort peaks in-place based on their center positions (x values).

        Parameters
        ----------
        order : str
            Sorting direction. Use:
                - "asc" for ascending (default)
                - "desc" for descending
        """
        reverse = {"asc": False, "desc": True}.get(order)
        if reverse is None:
            raise ValueError("order must be 'asc' or 'desc'")

        self._peaks.sort(key=lambda p: p['x'], reverse=reverse)
        return self

# ------------------------
# Generator
# ------------------------
class generator:
    def __init__(self, kind="gauss"):
        """define a peak generator gauss/lorentz/triangle"""
        self.kind = kind.lower()

    def __call__(self, x, x0, w, h):
        if self.kind.startswith("gauss"): #simplified Gaussian kernel
            # area: 0.6006*sqrt(pi)*w*h
            # standard deviation: 0.6006/sqrt(2)*w
            return h * np.exp(-((x - x0) / (0.6006 * w)) ** 2)
        elif self.kind.startswith("lorentz"):
            return  h * 1/(1+((x-x0) / (0.5*w)) ** 2)
        elif self.kind.startswith("triang"):
            d = np.abs(x - x0)
            return h * np.maximum(1 - d / (0.6006 * w), 0)
        else:
            raise ValueError(f"Unknown generator kind: {self.kind}")

    def __repr__(self):
        return f"<generator kind='{self.kind}'>"


# %% Example usage
if __name__ == "__main__":

    # output folder
    outputfolder = "./images" if os.path.isdir("./images") else ("../images" if os.path.isdir("../images") else None)

    # Define peaks
    p = peaks()
    p.add(x=10, w=2, h=1)
    p.add(x=20, w=2, h=1)
    p.add(x=30, w=2, h=1)

    # Convert to signal
    s = p.to_signal()  # defaults to x0=0, n=1000
    print(s)
    fig,_ = s.plot(); fig.print("simple_signal",outputfolder)

    # Generate variants
    s_noisy1 = s.add_noise(kind="gaussian", scale=0.01, bias=5)
    s_noisy2 = s.add_noise(kind="poisson", scale=10, bias="ramp")
    s_scaled = s * 0.5


    # Create a signal collection
    collection = signal_collection(s, s_noisy1, s_noisy2, s_scaled)
    print(collection)
    fig = collection.plot(); fig.print("collection_simple_signals",outputfolder)

    s_mean = collection.mean()
    s_sum = collection.sum()

    s_mean.plot(label="Mean")
    s_sum.plot(label="Sum")


    # Synthetic signals
    S,pS = signal_collection.generate_synthetic(
        n_signals=12,
        n_peaks=1,
        kinds=("gauss",),
        width_range=(0.5, 3),
        height_range=(1.0, 5.0),
        x_range=(0, 500),
        n_points=2048,
        normalize=False,
        seed=123
    )
    pS = peaks(pS) # flatten all peaks
    S.plot(fontsize=14)
    Sfull = S.mean()
    fig,_ = Sfull.plot(newfig=True); fig.print("synthetic_signal",outputfolder)

    # Transform
    dna = DNAsignal(Sfull)
    dna.compute_cwt()
    dna.transforms.plot().print("synthetic_transformed_signal",outputfolder)
    dna.encode_dna()
    dna.encode_dna_full()
    dna.plot_codes(4).print("synthetic_DNA_signal",outputfolder)
    dna.plot_scalogram().print("synthetic_scalogram",outputfolder)
    A=dna.codesfull[4]
    B=dna.codesfull[2]
    # Alignement
    A.align(B)
    A.html_alignment()
    print(A.wrapped_alignment(40))
    A.plot_alignment()
    As = A.to_signal()
    Bs = B.to_signal()
    ABs = signal_collection(As,Bs)
    plt.figure(); As.plot()
    ABs.plot()
    # substring (peaks) extraction
    pA_list = A.find("YAZB")
    pA =[s.to_signal() for s in pA_list]
    pAs = signal_collection(*[s.to_signal() for s in pA_list])
    pAs.plot()

    # Sinusoidal enoding
    dna.sinencode_dna()
    dna.sinencode_dna_full()
    dna.codes.plot().print("synthetic_DNA_raw",outputfolder)
    dna.sincodes.plot().print("synthetic_DNA_raw_sinencoded",outputfolder)
    dna.codesfull.plot().print("synthetic_DNA_letters",outputfolder)
    dna.sincodesfull.plot().print("synthetic_DNA_letters_sinencoded",outputfolder)

    # Signal mixtures and their classification
    Smix, pSmix, idSmix = signal_collection.generate_mixtures(
        n_mixtures=30,
        max_peaks=16,
        peaks_per_mixture=(4, 8),
        amplitude_range=(0.5, 2),
        n_signals=1,
        kinds=("gauss",),
        width_range=(0.5, 3),
        height_range=(1.0, 5.0),
        x_range=(0, 500),
        n_points=2048,
        normalize=False,
        seed=123
        )

    dnaSmix = Smix._toDNA(scales=[1,2,4,8,16,32])
    D = DNAsignal._pairwiseEntropyDistance(dnaSmix,scale=4,engine="bio")
    D.name = "Excess Entropy"
    Ddhalf, figD1 = D.dimension_variance_curve(); figD1.print("Entropy_dimensions",outputfolder)
    D.select_dimensions(10)
    figD2 = D.plot_dendrogram(); figD2.print("Entropy_dendrogram",outputfolder)
    figD3 = D.scatter3d(n_clusters=5); figD3.print("Entropy_scatter",outputfolder)

    J = DNAsignal._pairwiseJaccardMotifDistance(dnaSmix,scale=4)
    J.name = "YAZB Jaccard"
    Jdhalf,figJ1 = J.dimension_variance_curve(); figJ1.print("Jaccard_dimensions",outputfolder)
    J.select_dimensions(10)
    figJ2 = J.plot_dendrogram(); figJ2.print("Jaccard_dendrogram",outputfolder)
    figJ3 = J.scatter3d(n_clusters=5); figJ3.print("Jaccard_scatter3",outputfolder)

    L = DNAsignal._pairwiseLevenshteinDistance(dnaSmix,scale=4)
    L.name = "Levenshtein"
    Ldhalf,figL1 = L.dimension_variance_curve(); figL1.print("Lenvenshtein_dimensions",outputfolder)
    L.select_dimensions(10)
    figL2 = L.plot_dendrogram(); figL2.print("Lenvenshtein_dendrogram",outputfolder)
    figL3 = L.scatter3d(n_clusters=5); figL3.print("Lenvenshtein_scatter3",outputfolder)

    S = DNAsignal._pairwiseJensenShannonDistance(dnaSmix,scale=4)
    S.name = "Jensen-Shannon"
    Sdhalf,figJ1 = S.dimension_variance_curve(); figJ1.print("JensenShannon_dimensions",outputfolder)
    S.select_dimensions(10);
    figJ2 = S.plot_dendrogram(); figJ2.print("JensenShannon_dendrogram",outputfolder)
    figJ3 = S.scatter3d(n_clusters=5); figJ3.print("JensenShannon_scatter3",outputfolder)

    # Final summary table
    print('\nNumber of dimensions to recover 50% of total distance:\n')

    header = f"{'Distance Metric':<20} | {'dhalf':>6}"
    separator = "-" * len(header)
    rows = [
        ("Excess Entropy", Ddhalf),
        ("YAZB Jaccard", Jdhalf),
        ("Levenshtein", Ldhalf),
        ("Jensen-Shannon", Sdhalf),
    ]

    print(header)
    print(separator)
    for name, val in rows:
        print(f"{name:<20} | {val:6}")

    # Demo Alignment on Banner
    S1 = signal_collection.generate_synthetic(
        n_signals=8,
        n_peaks=2, kinds=("gauss",),
        width_range=(0.5, 2), height_range=(1.0, 5.0),
        x_range=(0, 127), n_points=128, normalize=False, seed=123)[0].mean()
    S1.name = "S1"

    S2 = signal_collection.generate_synthetic(
        n_signals=6,
        n_peaks=2, kinds=("gauss",),
        width_range=(0.4, 1.5), height_range=(1.0, 5.0),
        x_range=(0, 127), n_points=128, normalize=False, seed=456)[0].mean()
    S2.name = "S2"
    S12 = signal_collection(S1,S2,n=128)
    S12.name = "S1 vs. S2"
    S12.plot()
    S12dna = S12._toDNA(scales=[1,2,4,8,16,32])
    s1 = signal(y=S12dna[0].cwt_coeffs[2], name="S1")
    s2 =  signal(y=S12dna[1].cwt_coeffs[2], name="S2")
    s1code = S12dna[0].codesfull[2]
    s2code = S12dna[1].codesfull[2]
    s1code.align(s2code,"bio")
    S12dna_scale2 = signal_collection(s1,s2)
    S12dna_scale2.name="S1 vs. S2 - CWT at scale 2"
    S12dna_scale2.plot()
    s1code.plot_alignment()
    s1code.plot_mask()
    print(s1code.wrapped_alignment(width=1024))

    # --------------------- [SinEnCoder/DeCoder demo] --------------------------
    # 1. Construct a test input signal with smooth and jump segments
    x_smooth = np.linspace(0, 20, 100)
    x_jumps = np.array([25, 25, 26, 27, 100, 101, 130])
    x = np.concatenate([x_smooth, x_jumps])

    # 2. Initialize encoder with high d_model and N
    s = SinusoidalEncoder(d_model=128, N=10000)

    # 3. Fit auto-scaling to compress input into sinusoidal-friendly space
    s.fit_encoder(x, target_range=10)

    # 4. Encode and decode using all robust methods
    a = s.encode(x)
    decoded_lsq, err_lsq = s.decode(a, method='least_squares', return_error=True)
    decoded_svd, err_svd = s.decode(a, method='svd', return_error=True)

    # 5. Compare errors
    true = x
    lsq_error = np.abs(decoded_lsq - true)
    svd_error = np.abs(decoded_svd - true)

    # 6. Plot results
    fig, axs = plt.subplots(2, 2, figsize=(12, 8))
    axs[0, 0].plot(true, label="Original")
    axs[0, 0].plot(decoded_lsq, '--', label="Decoded (LSQ)")
    axs[0, 0].plot(decoded_svd, ':', label="Decoded (SVD)")
    axs[0, 0].set_title("Decoded vs Original")
    axs[0, 0].legend()

    axs[0, 1].plot(lsq_error, label="Abs Error (LSQ)")
    axs[0, 1].plot(svd_error, label="Abs Error (SVD)")
    axs[0, 1].set_yscale('log')
    axs[0, 1].set_title("Absolute Decoding Error (log scale)")
    axs[0, 1].legend()

    axs[1, 0].plot(err_lsq, label="Residual Norm (LSQ)")
    axs[1, 0].plot(err_svd, label="Residual Norm (SVD)")
    axs[1, 0].set_yscale('log')
    axs[1, 0].set_title("Reconstruction Residuals")
    axs[1, 0].legend()

    axs[1, 1].hist(lsq_error, bins=50, alpha=0.7, label="LSQ")
    axs[1, 1].hist(svd_error, bins=50, alpha=0.5, label="SVD")
    axs[1, 1].set_title("Histogram of Absolute Errors")
    axs[1, 1].legend()

    plt.suptitle("🌀 SinusoidalEncoder: Accuracy Evaluation", fontsize=14)
    plt.tight_layout()
    plt.show()
    fig.print("Sinusoidal-Encoder-Decoder",outputfolder)

    #%% naive Pseudoinverse example

    # Create a synthetic Gaussian signal
    x = np.linspace(0, 1, 1024)
    signal_data = np.exp(-(x - 0.5)**2 / (2 * 0.0025))

    # Choose a compatible continuous wavelet
    wavelet = pywt.ContinuousWavelet("mexh")  # Mexican Hat = Ricker

    # Scales for the transform
    scales = [1, 2, 4, 8, 16, 32]
    wavelet_matrix = []

    for scale in scales:
        cwt_result, _ = pywt.cwt(signal_data, [scale], wavelet)
        wavelet_matrix.append(cwt_result[0])

    # Stack transforms into matrix W
    W = np.vstack(wavelet_matrix)  # Shape: (n_scales, len(signal))

    # SVD-based pseudoinverse
    U, s, Vh = svd(W, full_matrices=False)
    rank = 6  # Truncate if desired
    S_inv = np.diag([1/s[i] if i < rank else 0 for i in range(len(s))])
    recon_signal = (U @ S_inv @ U.T @ W).sum(axis=0)

    # Plot
    figNaive = plt.figure(figsize=(8, 4))
    plt.plot(x, signal_data, label="Original", linewidth=2)
    plt.plot(x, recon_signal, '--', label="Reconstructed (pseudo-inv)", linewidth=2)
    plt.title("Pseudo-Inverse Reconstruction from Ricker CWT")
    plt.xlabel("x"); plt.ylabel("Amplitude"); plt.grid(True); plt.legend()
    plt.tight_layout(); plt.show()
    figNaive.print("Naive_PseudoInversion",outputfolder)


    # %% Advanced Pseudoinverse example

    # Synthetic signal with three Gaussians (including overlapping)
    x = np.linspace(0, 1, 1024)
    gauss1 = np.exp(-(x - 0.3)**2 / (2 * 0.002))
    gauss2 = np.exp(-(x - 0.5)**2 / (2 * 0.002))
    gauss3 = np.exp(-(x - 0.52)**2 / (2 * 0.002))  # Overlaps with gauss2
    signal_data = gauss1 + gauss2 + gauss3

    # Define scales and initialize
    scales = [1, 2, 4, 8, 16]
    wavelet_matrix = []
    energies = []

    # Use PyWavelets for energy normalization
    import pywt

    for scale in scales:
        coef, _ = pywt.cwt(signal_data, [scale], 'mexh')
        wavelet_matrix.append(coef[0])

        # Compute the wavelet kernel and its energy
        wavelet_kernel = pywt.ContinuousWavelet('mexh').wavefun(level=10)[0]
        energy = np.sqrt(np.sum(wavelet_kernel**2))  # L2 norm
        energies.append(energy)

    # Normalize each transformed scale by its wavelet energy
    W = np.vstack([w / e for w, e in zip(wavelet_matrix, energies)])

    # Perform SVD and pseudo-inverse reconstruction
    U, s, Vh = svd(W, full_matrices=False)
    S_inv = np.diag(1 / s)
    recon_signal = U.T @ W
    recon_signal = U @ S_inv @ recon_signal
    recon_signal = recon_signal.sum(axis=0)

    # Plot original and reconstructed signal
    figAdv = plt.figure(figsize=(10, 4))
    plt.plot(x, signal_data, label="Original Signal")
    plt.plot(x, recon_signal, label="Reconstructed (Weighted SVD)", linestyle="--")
    plt.xlabel("x")
    plt.ylabel("Amplitude")
    plt.title("Weighted Pseudoinverse Reconstruction from CWT")
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.show()
    figAdv.print("Weighted_PseudoInversion",outputfolder)

# %% Obsolete
"""
    # Early examples
    x = np.linspace(0, 1000, 1000)
    sig = np.exp(-((x - 500) ** 2) / (2 * 30 ** 2))

    s2d = sig2dna(sig)
    s2d.compute_cwt(scales=[1, 2, 4, 8])
    s2d.encode_dna(scale=4)

    print(s2d.get_code(scale=4))
    print("Entropy:", s2d.get_entropy(scale=4))

    matches = s2d.find_sequence("A", scale=4)
    print("Matches:", matches)

    reconstructed = s2d.reconstruct_signal(scale=4)

    plt.plot(sig, label='Original Signal')
    plt.plot(reconstructed, label='Reconstructed Signal', linestyle='--')
    plt.legend()
    plt.show()
"""
