The `align` package
===================

This package uses a pangenome as a reference to compute elements for a given genome, or a given set of proteins. As such, analysis that are usually run on multiple genomes can be run on the single genome or set of proteins that is provided. This subpackage depends on many of the other subpackages to run its analysis.
This package depends on the following packages:

- `formats`, to check the pangenome status.
- `annotate`, to read the given input files that can be gff or gbff.
- `cluster`, to write gene sequences from annotations.
- `RGP`, to eventually compute RGP and spot predictions.

It depends on the following modules:

- `pangenome`
- `utils`

Submodules
----------

ppanggolin.align.alignOnPang module
-----------------------------------

.. automodule:: ppanggolin.align.alignOnPang
   :members:
   :undoc-members:
   :show-inheritance: