===========================================================================
RegSeg: Structure-informed segmentation and registration of brain MR images
===========================================================================

.. image:: https://img.shields.io/badge/Citation-NeuroImage%20doi%3A10.1016%2Fj.neuroimage.2016.05.011-blue.svg
  :target: https://doi.org/10.1016/j.neuroimage.2016.05.011

.. image:: https://img.shields.io/badge/License-MIT-blue.svg
  :target: https://github.com/oesteban/RegSeg/blob/de89dfb01abed3778e9764ab12fdcdb2dfc187eb/LICENSE

RegSeg is a simultaneous segmentation and registration method that uses
active contours without edges (ACWE) extracted from structural images.
The contours evolve through a free-form deformation field supported by the
B-spline basis to optimally map the contours onto the data in the target
space.

.. image :: docs/static/graphical-abstract.png


.. topic:: **When using this software in your research, please credit the authors referencing the following paper:**

    Esteban O, Zosso D, Daducci A, Bach-Cuadra M, Ledesma-Carbayo MJ, Thiran JP, Santos A;
    *Surface-driven registration method for the structure-informed segmentation of diffusion MR images*;
    NeuroImage 139:450-461; 1 October 2016;
    doi:`10.1016/j.neuroimage.2016.05.011 <https://doi.org/10.1016/j.neuroimage.2016.05.011>`_.


----------------------
Experimental framework
----------------------

RegSeg is distributed along with the software instrumentation to benchmark it.
The experimental framework is written in Python and uses nipype.

We tested the functionality of regseg using four digital phantoms warped with
known and randomly generated deformations, where subvoxel accuracy was achieved.
We then applied regseg to a registration/segmentation task using 16 real diffusion MRI
datasets from the Human Connectome Project, which were warped by realistic and nonlinear
distortions that are typically present in these data.
We computed the misregistration error of the contours estimated by regseg with respect to
their theoretical location using the ground truth, thereby obtaining a 95% CI of 0.56–0.66
mm distance between corresponding mesh vertices, which was below the 1.25 mm isotropic
resolution of the images.
We also compared the performance of our proposed method with a widely used registration tool,
which showed that regseg outperformed this method in our settings.


------------
Installation
------------
::

  mkdir Release
  cd Release
  ccmake ../Code/ -G"Eclipse CDT4 - Unix Makefiles" -DCMAKE_BUILD_TYPE=Release -DITK_DIR=/usr/local/lib/cmake/ITK-4.7/


-----------
MIT License
-----------

.. include LICENSE
