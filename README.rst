:author: Lukas Turcani
:Documentation: https://stk.readthedocs.io
:Discord: https://discord.gg/zbCUzuxe2B

.. figure:: docs/source/figures/stk.png

.. image:: https://github.com/lukasturcani/stk/workflows/tests/badge.svg?branch=master
    :target: https://github.com/lukasturcani/stk/actions?query=branch%3Amaster

.. image:: https://readthedocs.org/projects/stk/badge/?version=latest
    :target: https://stk.readthedocs.io

Overview
========

``stk`` is a Python library which allows construction and
manipulation of complex molecules, as well as automatic
molecular design, and the creation of molecular, and molecular
property, databases. The documentation of ``stk`` is available on
https://stk.readthedocs.io and the project's Discord server can be
joined through https://discord.gg/YvwdcjKf.

Installation
============

To get ``stk``, you can install it with pip::

    $ pip install stk

Make sure you also install rdkit, which is a dependency::

    $ conda install -c conda-forge rdkit

If you would like to get updated when a new release of ``stk`` comes
out, which happens pretty regularly, click on the ``watch`` button on
the top right corner of the GitHub page. Then select ``Releases only``
from the dropdown menu.

You can see the latest releases here:

    https://github.com/lukasturcani/stk/releases

There will be a corresponding release on ``pip`` for each release
on GitHub, and you can update your ``stk`` with::

    $ pip install stk --upgrade

How To Cite
===========

If you use ``stk`` please cite

    https://github.com/lukasturcani/stk

and

    https://aip.scitation.org/doi/10.1063/5.0049708


Publications
============

about stk
---------

* `stk: An Extendable Python Framework for Automated Molecular and
  Supramolecular Structure Assembly and Discovery`__

__ https://aip.scitation.org/doi/10.1063/5.0049708

* (Out of date) `stk: A Python Toolkit for Supramolecular Assembly`__
  | chemrxiv__

__ https://onlinelibrary.wiley.com/doi/abs/10.1002/jcc.25377
__ https://chemrxiv.org/articles/STK_A_Python_Toolkit_for_Supramolecular_Assembly/6127826

using stk
---------

* `High-throughput Computational Evaluation of Low Symmetry Pd2L4
  Cages to Aid in System Design`__

__ https://onlinelibrary.wiley.com/doi/10.1002/anie.202106721

* `Forecasting System of Computational Time of DFT/TDDFT Calculations
  under the Multiverse Ansatz via Machine Learning and
  Cheminformatics`__

__ https://pubs.acs.org/doi/full/10.1021/acsomega.0c04981

* `Using High-throughput Virtual Screening to Explore the
  Optoelectronic Property Space of Organic Dyes; Finding
  Diketopyrrolopyrrole Dyes for Dye-sensitized Water Splitting and
  Solar Cells`__

__ https://pubs.rsc.org/en/content/articlelanding/2021/SE/D0SE00985G#!divAbstract

* `Accelerated Discovery of Organic Polymer Photocatalysts for Hydrogen
  Evolution from Water through the Integration of Experiment and
  Theory`__

__ https://pubs.acs.org/doi/abs/10.1021/jacs.9b03591

* `Structurally Diverse Covalent Triazine-Based Framework Materials for
  Photocatalytic Hydrogen Evolution from Water`__

__ https://pubs.acs.org/doi/full/10.1021/acs.chemmater.9b02825

* `Mapping Binary Copolymer Property Space with Neural Networks`__

__ https://pubs.rsc.org/ko/content/articlehtml/2019/sc/c8sc05710a

* `An Evolutionary Algorithm for the Discovery of Porous Organic
  Cages`__ | chemrxiv__

__ https://pubs.rsc.org/en/content/articlelanding/2018/sc/c8sc03560a#!divAbstract
__ https://chemrxiv.org/articles/An_Evolutionary_Algorithm_for_the_Discovery_of_Porous_Organic_Cages/6954557

* `Machine Learning for Organic Cage Property Prediction`__
  | chemrxiv__

__ https://pubs.acs.org/doi/10.1021/acs.chemmater.8b03572
__ https://chemrxiv.org/articles/Machine_Learning_for_Organic_Cage_Property_Prediction/6995018

* `A High-Throughput Screening Approach for the Optoelectronic
  Properties of Conjugated Polymers`__ | chemrxiv__

__ https://pubs.acs.org/doi/abs/10.1021/acs.jcim.8b00256
__ https://chemrxiv.org/articles/A_High-Throughput_Screening_Approach_for_the_Optoelectronic_Properties_of_Conjugated_Polymers/6181841

* `Computationally-Inspired Discovery of an Unsymmetrical Porous
  Organic Cage`__ | chemrxiv__

__ https://pubs.rsc.org/en/content/articlelanding/2018/nr/c8nr06868b#!divAbstract
__ https://chemrxiv.org/articles/Computationally-Inspired_Discovery_of_an_Unsymmetrical_Porous_Organic_Cage/6863684

* `Maximising the Hydrogen Evolution Activity in Organic Photocatalysts
  by co-Polymerisation`__

__ https://pubs.rsc.org/en/Content/ArticleLanding/TA/2018/C8TA04186E#!divAbstract


Acknowledgements
================

I began developing this code when I was working in the Jelfs group,
http://www.jelfs-group.org/, whose members often provide me with
very valuable feedback, which I gratefully acknowledge.
