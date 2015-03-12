=================
ACWE-Registration
=================

This project is an extension of the Active Contours Without Edges (ACWE) framework to preform
non-linear registration of strong priors to a reference image with subvoxel precision.


-------
Credits
-------

Authors
-------

* Oscar Esteban (@oesteban): registration and evaluation methods implementation,
  design of evaluation methods.
* Dominique Zosso (@zosso): registration method formulation
* Alessandro Daducci (@daducci): early datasets generation, interpretation of results.
* Meritxell Bach-Cuadra (@meribach): interpretation of results, project coordination.

Supervision
-----------

* Andrés Santos Lleó
* Jean-Philippe Thiran
* M.-J. Ledesma-Carbayo


-----------
Compilation
-----------
::

  mkdir Debug
  cd Debug
  ccmake ../Code/ -G"Eclipse CDT4 - Unix Makefiles" -DCMAKE_BUILD_TYPE=Debug -DITK_DIR=/usr/local/lib/cmake/ITK-4.2/

