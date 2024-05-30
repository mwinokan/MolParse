.. MolParse documentation master file, created by
   sphinx-quickstart on Wed Mar 13 13:02:45 2024.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

======================
MolParse Documentation
======================

**Warning: this documentation is barebones**

MolParse is a python package for parsing, modifying, and analysis of molecular structure files built on top of ASE.

Installation
============

Install from PyPI

::

   pip install --upgrade molparse

Getting Started
===============

To parse a file try:

::

import molparse as mp
sys = mp.parse(file)

This will return a hierarchical :class:`.System` object (for .pdb and .gro files) which is structured as follows:

::

   System
   |__ Chain
       |__ Residue
           |__ Atom

N.B. All groups of Atom objects inherit from AtomGroup

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   API Reference <api_reference>
   .. molparse
   .. io
   .. group
   .. system
   .. chain
   .. residue
   .. atom




.. Indices and tables
.. ==================

.. * :ref:`genindex`
.. * :ref:`modindex`
.. * :ref:`search`
