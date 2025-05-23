#!/usr/bin/env python3
# -*- coding: utf-8 -*-


__all__ = ['PrintableFigure', 'custom_plt_figure', 'custom_plt_subplots', 'is_valid_figure', 'print_figure', 'print_pdf', 'print_png', 'print_svg']
"""
📁 figprint.py — Utilities for saving and displaying Matplotlib figures with printing options.

This module extends Matplotlib's `Figure` class by adding convenient methods to export
figures in PDF, PNG, and SVG formats. It also overrides `plt.figure` and `plt.subplots`
to return enhanced `PrintableFigure` objects by default.

Typical usage:
--------------
    import matplotlib.pyplot as plt
    from figprint import print_figure

    # Automatically uses PrintableFigure
    fig, ax = plt.subplots()
    ax.plot([0, 1], [0, 1])
    fig.print()  # Saves to PNG and PDF

    # Advanced control
    fig.print_png("myplot", overwrite=True)
    fig.print_svg("myplot", overwrite=True)

Example:
--------
# myplotmodule.py
from matplotlib import pyplot as plt
import figprint  # Automatically patches pyplot

def example_plot():
    fig, ax = plt.subplots()
    ax.plot([1, 2, 3], [4, 5, 6])
    ax.set_title("Example Plot")
    fig.print("example", overwrite=True)


Author: Generative Simulation Initiative/Olivier Vitrac, PhD
Created: 2025-05-12
"""

import os
from datetime import datetime
import matplotlib.pyplot as plt
from matplotlib.figure import Figure

_fig_metadata_atrr_ = "filename"


def is_valid_figure(fig):
    """
    Check if `fig` is a valid and open Matplotlib figure.
    Parameters:
        fig (object): Any object to check.

    Returns:
        bool: True if the object is an open Matplotlib Figure.
    """
    return isinstance(fig, Figure) and hasattr(fig, 'canvas') and fig.canvas is not None


def _generate_figname(fig, extension):
    """
    Generate a cleaned filename using figure metadata or the current datetime.

    Parameters:
        fig (Figure): Matplotlib figure object.
        extension (str): File extension (e.g. '.pdf').

    Returns:
        str: Cleaned filename with correct extension.
    """
    if hasattr(fig, _fig_metadata_atrr_):
        filename = getattr(fig, _fig_metadata_atrr_)
    else:
        filename = "fig" + datetime.now().strftime("%Y%m%d_%H%M%S")
    filename = filename.strip().replace(" ", "_")
    if not filename.lower().endswith(extension):
        filename += extension
    return filename


def _print_generic(fig, extension, filename, destinationfolder, overwrite, dpi):
    """
    Generic saving logic for figure files.

    Parameters:
        fig (Figure): The figure to save.
        extension (str): File extension (e.g. '.pdf').
        filename (str): Optional filename.
        destinationfolder (str): Folder path.
        overwrite (bool): Overwrite existing file?
        dpi (int): Resolution (ignored for SVG).
    """
    if not is_valid_figure(fig):
        print("no valid figure")
        return

    filename = filename or _generate_figname(fig, extension)
    if not filename.endswith(extension):
        filename += extension
    filepath = os.path.join(destinationfolder, filename)

    if not overwrite and os.path.exists(filepath):
        print(f"File {filepath} already exists. Use overwrite=True to replace it.")
        return

    fig.savefig(filepath, format=extension.lstrip("."), dpi=None if extension == ".svg" else dpi, bbox_inches="tight")
    print(f"Saved {extension.upper()[1:]}: {filepath}")


def print_pdf(fig, filename="", destinationfolder=os.getcwd(), overwrite=False, dpi=300):
    """
    Save a figure as a PDF.
    Example:
    --------
    >>> fig, ax = plt.subplots()
    >>> ax.plot([0, 1], [0, 1])
    >>> print_pdf(fig, "myplot", overwrite=True)
    """
    _print_generic(fig, ".pdf", filename, destinationfolder, overwrite, dpi)


def print_png(fig, filename="", destinationfolder=os.getcwd(), overwrite=False, dpi=300):
    """
    Save a figure as a PNG.

    Example:
    --------
    >>> print_png(fig, "myplot", overwrite=True)
    """
    _print_generic(fig, ".png", filename, destinationfolder, overwrite, dpi)


def print_svg(fig, filename="", destinationfolder=os.getcwd(), overwrite=False):
    """
    Save a figure as an SVG (vector format, dpi-independent).

    Example:
    --------
    >>> print_svg(fig, "myplot", overwrite=True)
    """
    _print_generic(fig, ".svg", filename, destinationfolder, overwrite, dpi=None)


def print_figure(fig, filename="", destinationfolder=os.getcwd(), overwrite=False,
 dpi={"png": 150, "pdf": 300, "svg": None}, what = ["pdf","png","svg"]):
    """
    Save a figure as PDF, PNG, and SVG.

    Example:
    --------
    >>> print_figure(fig, "full_output", overwrite=True)

    Parameters:
        fig (Figure): Figure to save.
        filename (str): Optional base filename.
        destinationfolder (str): Folder path.
        overwrite (bool): Whether to overwrite existing files.
        dpi (dict): Dictionary of resolution per format.
        what (list): list what to print, default = ["pdf","png","svg"]
    """
    if not isinstance(what,(list,tuple)):
        what = [what]
    if is_valid_figure(fig):
        if "pdf" in what: print_pdf(fig, filename, destinationfolder, overwrite, dpi["pdf"])
        if "png" in what: print_png(fig, filename, destinationfolder, overwrite, dpi["png"])
        if "svg" in what: print_svg(fig, filename, destinationfolder, overwrite)
    else:
        print("no valid figure")


class PrintableFigure(Figure):
    """
    Enhanced Matplotlib Figure with custom show and export methods.

    Example:
    --------
    >>> fig, ax = plt.subplots()
    >>> ax.plot([0, 1], [1, 0])
    >>> fig.print("diag1")  # Saves PDF, PNG, SVG
    """

    def show(self, display_mode=None):
        """
        Display figure intelligently based on context (Jupyter/script).

        Parameters:
            display_mode (str): 'auto' or 'classic' (default is 'auto').
        """
        try:
            get_ipython
            if display_mode is None or display_mode == "auto":
                display(self)
            else:
                super().show()
        except NameError:
            super().show()

    def print(self, filename="", destinationfolder=os.getcwd(), overwrite=True,
            dpi={"png": 150, "pdf": 300, "svg": None}):
        """
        Save figure in PDF, PNG, and SVG formats.

        Example:
        --------
        >>> fig.print("summary_figure")
        """
        print_figure(self, filename, destinationfolder, overwrite, dpi)

    def print_pdf(self, filename="", destinationfolder=os.getcwd(), overwrite=False, dpi=300):
        print_pdf(self, filename, destinationfolder, overwrite, dpi)

    def print_png(self, filename="", destinationfolder=os.getcwd(), overwrite=False, dpi=300):
        print_png(self, filename, destinationfolder, overwrite, dpi)

    def print_svg(self, filename="", destinationfolder=os.getcwd(), overwrite=False):
        print_svg(self, filename, destinationfolder, overwrite)


# Save original constructor references
original_plt_figure = plt.figure
original_plt_subplots = plt.subplots

def custom_plt_figure(*args, **kwargs):
    """
    Override `plt.figure()` to return PrintableFigure by default.

    Returns:
        PrintableFigure
    """
    kwargs.setdefault("FigureClass", PrintableFigure)
    return original_plt_figure(*args, **kwargs)

def custom_plt_subplots(*args, **kwargs):
    """
    Override `plt.subplots()` to return PrintableFigure.

    Returns:
        (PrintableFigure, Axes)
    """
    fig, ax = original_plt_subplots(*args, **kwargs)
    fig.__class__ = PrintableFigure
    return fig, ax

# Apply overrides globally
plt.figure = custom_plt_figure
plt.subplots = custom_plt_subplots
plt.FigureClass = PrintableFigure
plt.rcParams['figure.figsize'] = (8, 6)
