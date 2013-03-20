=================
ACWE-Registration
=================

This project is an extension of the Active Contours Without Edges (ACWE) framework to preform
non-linear registration of strong priors to a reference image with subvoxel precision.


-------
Credits
-------

* Oscar Esteban (oesteban): implementation
* Dominique Zosso (zosso): mathematical formulation
* Jean-Philippe Thiran
* Meritxell Bach-Cuadra
* Andrés Santos Lleó



-----------
Compilation
-----------
::

  mkdir Debug
  cd Debug
  ccmake ../Code/ -G"Eclipse CDT4 - Unix Makefiles" -DCMAKE_BUILD_TYPE=Debug -DITK_DIR=/usr/local/lib/cmake/ITK-4.2/

