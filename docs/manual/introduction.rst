.. _introduction:

Introduction
============

This is the documentation for the implementation of a single-shot emittance measurement at FACET.

Prerequisites
-------------

Python 3
^^^^^^^^

:mod:`oneshot` works with Python 3 and up, which should be installed via apt-get on \*nix, `Macports <https://www.macports.org/>`_ on Apple machines, or downloaded from https://www.python.org/downloads/.

NumPy and SciPy
^^^^^^^^^^^^^^^

:mod:`oneshot` depends on `NumPy <http://www.numpy.org/>`_ and `SciPy <http://www.scipy.org/>`_ to manipulate data.

`NumPy <http://www.numpy.org/>`_ has dependencies such as BLAS, LAPACK, and ATLAS, which makes downloading building form source is difficult. Installation via apt-get or `Macports <https://www.macports.org/>`_ is highly recommended in order to handle these dependencies. It is possible to `download <http://www.scipy.org/scipylib/download.html>`_ or to `build from source <http://www.scipy.org/scipylib/building/index.html#building>`_.

`SciPy <http://www.scipy.org/>`_ is similar to `NumPy <http://www.numpy.org/>`_, and it is usually easiest to install the two at the same time - they have nearly identical install methods, to the point that at times it is difficult to tell the two apart. It is easiest to follow the `install instructions <http://www.scipy.org/install.html>`_.
