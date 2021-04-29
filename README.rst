:author: Lukas Turcani
:Documentation: https://stk.readthedocs.io

.. figure:: docs/source/figures/stk.png

.. image:: https://github.com/lukasturcani/stk/workflows/tests/badge.svg?branch=master
    :target: https://github.com/lukasturcani/stk/actions?query=branch%3Amaster

.. image:: https://readthedocs.org/projects/stk/badge/?version=latest
    :target: https://stk.readthedocs.io

Overview
========

``stk`` is a Python library which allows construction and
manipulation of complex molecules, as well as automatic
molecular design, and the creation of molecular, and molecular property,
databases.

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

    https://onlinelibrary.wiley.com/doi/abs/10.1002/jcc.25377


Publications
============

about stk
---------

* `stk: A Python Toolkit for Supramolecular Assembly`_ | chemrxiv__

__ https://chemrxiv.org/articles/STK_A_Python_Toolkit_for_Supramolecular_Assembly/6127826

.. _`stk: A Python Toolkit for Supramolecular Assembly`: https://onlinelibrary.wiley.com/doi/abs/10.1002/jcc.25377

using stk
---------

* `An Evolutionary Algorithm for the Discovery of Porous Organic Cages`_ | chemrxiv__

__ https://chemrxiv.org/articles/An_Evolutionary_Algorithm_for_the_Discovery_of_Porous_Organic_Cages/6954557
.. _`An Evolutionary Algorithm for the Discovery of Porous Organic Cages`: https://pubs.rsc.org/en/content/articlelanding/2018/sc/c8sc03560a#!divAbstract

* `Machine Learning for Organic Cage Property Prediction`_ | chemrxiv__

__ https://chemrxiv.org/articles/Machine_Learning_for_Organic_Cage_Property_Prediction/6995018
.. _`Machine Learning for Organic Cage Property Prediction`: https://pubs.acs.org/doi/10.1021/acs.chemmater.8b03572


* `A High-Throughput Screening Approach for the Optoelectronic Properties of Conjugated Polymers`_ | chemrxiv__

__ https://chemrxiv.org/articles/A_High-Throughput_Screening_Approach_for_the_Optoelectronic_Properties_of_Conjugated_Polymers/6181841
.. _`A High-Throughput Screening Approach for the Optoelectronic Properties of Conjugated Polymers`: https://pubs.acs.org/doi/abs/10.1021/acs.jcim.8b00256

* `Computationally-Inspired Discovery of an Unsymmetrical Porous Organic Cage`_ | chemrxiv__

__ https://chemrxiv.org/articles/Computationally-Inspired_Discovery_of_an_Unsymmetrical_Porous_Organic_Cage/6863684
.. _`Computationally-Inspired Discovery of an Unsymmetrical Porous Organic Cage`: https://pubs.rsc.org/en/content/articlelanding/2018/nr/c8nr06868b#!divAbstract

* `Maximising the Hydrogen Evolution Activity in Organic Photocatalysts by co-Polymerisation`_

.. _`Maximising the Hydrogen Evolution Activity in Organic Photocatalysts by co-Polymerisation`: https://pubs.rsc.org/en/Content/ArticleLanding/TA/2018/C8TA04186E#!divAbstract


Acknowledgements
================

I began developing this code when I was working in the Jelfs group,
http://www.jelfs-group.org/, whose members often provide me with
very valuable feedback, which I gratefully acknowledge.
