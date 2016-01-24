#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: oesteban
# @Date:   2015-11-19 16:44:27
# @Last Modified by:   oesteban
# @Last Modified time: 2016-01-18 20:50:16


def main():
    from glob import glob
    from setuptools import setup

    setup(
        name='pyacwereg',
        version=0.1,
        description='',
        author_email='code@oscaresteban.es',
        url='https://github.com/oesteban/RegSeg',
        download_url='',
        license='GPL',
        packages=['pyacwereg', 'pyacwereg.data', 'pyacwereg.interfaces',
                  'pyacwereg.workflows'],
        package_data={'pyacwereg': ['data/*.json', 'data/*.txt']},
        scripts=['tools/extract_hcp.py', 'tools/run_evaluations.py',
                 'tools/run_phantoms.py'],
        install_requires=["nipype", "nibabel", "pandas", "seaborn"],
        zip_safe=False)

if __name__ == "__main__":
    import os
    import sys

    local_path = os.path.dirname(os.path.abspath(sys.argv[0]))
    os.chdir(local_path)
    sys.path.insert(0, local_path)

    main()
