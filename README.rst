===========================================================================
RegSeg: Structure-informed segmentation and registration of brain MR images
===========================================================================


RegSeg is a simultaneous segmentation and registration method that uses
active contours without edges (ACWE) extracted from structural images.
The contours evolve through a free-form deformation field supported by the
B-spline basis to optimally map the contours onto the data in the target
space.

.. image:: https://zenodo.org/badge/doi/10.5281/zenodo.17348.svg   :target: http://dx.doi.org/10.5281/zenodo.17348


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


-------
Credits
-------


Software
--------

https://github.com/oesteban/RegSeg/graphs/contributors


Contributions
-------------

* Oscar Esteban (@oesteban)
* Dominique Zosso (@zosso)
* Alessandro Daducci (@daducci)
* Meritxell Bach-Cuadra (@meribach)
* Andrés Santos Lleó (@asantos5)
* Jean-Philippe Thiran (@thiran)
* M.-J. Ledesma-Carbayo (@mledesmacarbayo)

All the authors contributed to this study.
OE implemented the method, designed and conducted the experiments, wrote the paper,
simulated the phantoms, and prepared the real data.
DZ devised and drafted the registration method, generated early phantom datasets,
and collaborated in the implementation of the method.
AD, MBC, and MJLC interpreted the results.
AD, MBC, MJLC, JPT, and AS advised on all aspects of the study.


-------
License
-------

RegSeg is free software: you can redistribute it and/or modify it under the terms of the
`GNU General Public License <http://www.gnu.org/copyleft/gpl.html>`_ as published by the
`Free Software Foundation <http://www.fsf.org/>`_, either version 3 of the License, or
(at your option) any later version.

RegSeg is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
See the `GNU General Public License <http://www.gnu.org/copyleft/gpl.html>`_ for more details.
You should have received a copy of `GNU General Public License <http://www.gnu.org/copyleft/gpl.html>`_
along with RegSeg. If not, see http://www.gnu.org/licenses/.
