Running the experiments
=======================

A brief guide to run the experiments


Real Data
---------

.. highlight:: bash

   run_evaluations.py -N NeuroImagePipeline -S ../ -s $( cat subject_selection.txt ) --nthreads 8 --log_dir mylogs &> mylogs/ni2015.001
