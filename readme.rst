:author: Lukas Turcani

.. figure:: docs/source/figures/stk.png

Requires Python 3.6.

This code was developed while doing my PhD in the Jelfs group,
http://www.jelfs-group.org/.

Overview
========

``stk`` is a Python library for building, manipulating, optimizing and
designing molecules, including a genetic algorithm for the
automated design of molecules. For more information, the documentation can be found at
https://lukasturcani.github.io/stk/docs/build/html.

Installation
============

There are two ways to install ``stk``. The first is easier, using
pip. The second involves downloading the source code form GitHub. Using
GitHub is the best way to get the most up-to-date code.

In both cases ``rdkit`` needs to be installed, see `installing rdkit`_.

pip
---

In the terminal run

.. code-block:: bash

    pip install stk

GitHub
------

Installation is simple and has two parts. First, installing ``stk``, second,
installing any libraries it requires.

installing stk
--------------

1. Download ``stk`` from GitHub by clicking on the green "Clone or download"
   button and then on "Download ZIP".
2. Extract the downloaded zip file.
3. Using the terminal, go into the extracted directory and run:

.. code-block:: bash

    python setup.py install

4. ``stk`` should be installed, you can delete the downloaded and
   extracted folders.

installing dependencies
-----------------------

Apart from ``rdkit`` all dependcies of ``stk`` are listed in
``requirements.txt``. ``rdkit`` cannot be installed via pip so it is
listed separately.

rdkit
.....

See `installing rdkit`_.

others
......

All other requirements can be installed with pip.

.. code-block:: bash

    pip install -r requirements.txt

installing rdkit
----------------

The simplest way to get ``rdkit`` is to use the Anaconda distribution of
Python, which can be found on https://www.anaconda.com/download/. It
comes with a lot of scientific libraries already installed and the
conda package manager.

After installing Anaconda run

.. code-block:: bash

    conda install -c rdkit rdkit

If you do not want to use Anaconda, the best place to look for advice
on installing ``rdkit`` is https://github.com/rdkit/rdkit.

Publications
============

about stk
---------

* `stk: A Python Toolkit for Supramolecular Assembly`_ | `chemrxiv`_

.. _`stk: A Python Toolkit for Supramolecular Assembly`: https://onlinelibrary.wiley.com/doi/abs/10.1002/jcc.25377
.. _`chemrxiv`: https://chemrxiv.org/articles/STK_A_Python_Toolkit_for_Supramolecular_Assembly/6127826

using stk
---------

* `An Evolutionary Algorithm for the Discovery of Porous Organic Cages`_ | `chemrxiv`_

.. _`An Evolutionary Algorithm for the Discovery of Porous Organic Cages`: https://pubs.rsc.org/en/content/articlelanding/2018/sc/c8sc03560a#!divAbstract
.. _`chemrxiv`: https://chemrxiv.org/articles/An_Evolutionary_Algorithm_for_the_Discovery_of_Porous_Organic_Cages/6954557

* `Machine Learning for Organic Cage Property Prediction`_ | `chemrxiv`_

.. _`Machine Learning for Organic Cage Property Prediction`:
.. _`chemrxiv`: https://chemrxiv.org/articles/Machine_Learning_for_Organic_Cage_Property_Prediction/6995018

* `A High-Throughput Screening Approach for the Optoelectronic Properties of Conjugated Polymers`_ | `chemrxiv`_

.. _`A High-Throughput Screening Approach for the Optoelectronic Properties of Conjugated Polymers`: https://pubs.acs.org/doi/abs/10.1021/acs.jcim.8b00256
.. _`chemrxvi`: https://chemrxiv.org/articles/A_High-Throughput_Screening_Approach_for_the_Optoelectronic_Properties_of_Conjugated_Polymers/6181841

* `Computationally-Inspired Discovery of an Unsymmetrical Porous Organic Cage`_ | `chemrxiv`_

.. _`Computationally-Inspired Discovery of an Unsymmetrical Porous Organic Cage`:
.. _`chemrxiv`: https://chemrxiv.org/articles/Computationally-Inspired_Discovery_of_an_Unsymmetrical_Porous_Organic_Cage/6863684

* `Maximising the Hydrogen Evolution Activity in Organic Photocatalysts by co-Polymerisation`_

.. _`Maximising the Hydrogen Evolution Activity in Organic Photocatalysts by co-Polymerisation`: https://pubs.rsc.org/en/Content/ArticleLanding/TA/2018/C8TA04186E#!divAbstract
