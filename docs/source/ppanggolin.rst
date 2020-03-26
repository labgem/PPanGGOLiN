Developper doc
==============

``ppanggolin`` is both a command line tool and a python library for comparative genomics. It tries to prodive a solution for using cutting-edge methods for large scale comparative analysis and stores any computed results in a compact format so that they can be reused at will.
This part of the documentation is made for people that want to use PPanGGOLiN as a python library, or for those that need to maintain the package or want to modify it.

If you were looking for the command line tool documentation of PPanGGOLiN, you should check the github wiki instead.

Subpackages
-----------

There is a ppanggolin subpackage for each specific step of the analysis. Each subpackage is associated to one or more subcommand.

.. toctree::
   :maxdepth: 1

   ppanggolin.RGP
   ppanggolin.align
   ppanggolin.annotate
   ppanggolin.cluster
   ppanggolin.figures
   ppanggolin.formats
   ppanggolin.graph
   ppanggolin.info
   ppanggolin.nem
   ppanggolin.workflow

Submodules
----------

Submodules includes all of the basic classes of PPanGGOLiN that will be used by the subpackages.

.. toctree::
   :maxdepth: 2

   classes/pangenome
   classes/edge
   classes/geneFamily

ppanggolin.genome module
------------------------

.. automodule:: ppanggolin.genome
   :members:
   :undoc-members:
   :show-inheritance:

ppanggolin.main module
----------------------

.. automodule:: ppanggolin.main
   :members:
   :undoc-members:
   :show-inheritance:

ppanggolin.region module
------------------------

.. automodule:: ppanggolin.region
   :members:
   :undoc-members:
   :show-inheritance:

ppanggolin.utils module
-----------------------

.. automodule:: ppanggolin.utils
   :members:
   :undoc-members:
   :show-inheritance:


Module contents
---------------

.. automodule:: ppanggolin
   :members:
   :undoc-members:
   :show-inheritance:
