.. ParGeo documentation master file, created by
   sphinx-quickstart on Mon Feb 12 16:24:42 2024.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to ParGeo's documentation!
==================================

`Overview <#overview>`__ \| `Installation <#installation>`__ \| `Basic
Usage <#basic-usage>`__ \| `Cite <#cite>`__ \| `Contents <#contents>`__ \| `Indices and tables <#indices-and-tables>`__

Generate complex two-dimensional geometries for simulations based on
mesh discretizations.

Overview
========

With ParGeo you can create complex *domains*, where a domain is a
two-dimensional geometry which is described by a collection of shapely
`Polygons <https://shapely.readthedocs.io/en/stable/reference/shapely.Polygon.html#shapely.Polygon>`__
and
`MultiPolygons <https://shapely.readthedocs.io/en/stable/reference/shapely.MultiPolygon.html#shapely.MultiPolygon>`__
called *subdomains*.

This domain could then be used for mesh generation, to further perform
2D mesh based simulations, like the Finite Element Methods (FEM). We
provide mesh generating functionalities based on
`Gmsh <https://gmsh.info>`__.

To generate the desired domain you will use the ``pargeo.domain.Domain``
class and sequentially add subdomains. Each subdomain is associated with
a visibility level. The visibility level is an integer, and subdomains
with a higher visibility level will cover those with a lower one and be
merged with those having the same one.

Following features support you with the domain-generating process:

-  **Create Geometries** (``pargeo.geometry``):

   Create common geometries that can easly be transformed into shapely
   Polygons.

   Following geometries are provided:

   -  Rectangle
   -  Circle / Ellipse
   -  Stellar
   -  NStar
   -  Raindrop

-  **Transform Geometries** (``pargeo.transform``)

   Pass transforms to your domain while adding new sub-domains. The
   sub-domain will then be first accordingly transformed before being
   added to the domain.

   Following transforms are provided:

   -  Repeat
   -  Periodic
   -  Diffeomorphism

-  **Impose Constraints** (``pargeo.constraint``)

   Pass constraints to the domain while adding new sub-domains. The new
   sub-domain is only added, if the constraints are met.

   Following constraints are provided:

   -  DistanceConstraint

Check out the
`documentation <https://pargeo.readthedocs.io/en/latest/>`__ for more
information.

Installation
------------

Install the latest version of ParGeo from PyPI using pip:

.. code:: bash

   pip install pargeo

Basic Usage
-----------

Here’s a quick example of how to use ParGeo.

Cite
----

If you use ParGeo in your research, please cite it. You can use the
following BibTeX entry:

.. code:: bibtex

   @software{Schafer_ParGeo_2024,
   author = {Schäfer, Till and Gruhlke, Robert},
   month = feb,
   title = {{ParGeo}},
   url = {https://github.com/Priusds/ParGeo},
   version = {0.2.0},
   year = {2024}
   }

Contents
--------

.. toctree::
   :maxdepth: 1

   modules
   notebooks/getting-started

Indices and tables
------------------

* :ref:`genindex`
* :ref:`modindex`

