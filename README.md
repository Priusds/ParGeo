# Bubbles

With `Bubbles` you can simulate random porose materials as GMSH objects, that can be further used in FEniCS.
The porose material can be two-dimensional or three-dimensional. _TODO: Add some images to clarify what is possible with
this software._

## Installation

Install the package by running
```bash
(pyvenv) $ pip install -e .
```
in a virtual python environment.

Tested with:
- macOS Ventura Version 13.3.1 with Apple M2 Chip
- Gmsh 4.11.1
- Dolfin 2019.1.0
- Python 3.11.1

## Documentation

### Create Dolfin readable mesh
Once the `.GEO` file is written, following commands can be used to
create a dolfin compatible `.XML` file representing the mesh.

Write the mesh as a `.MSH` file, in the case of a 2d mesh
```bash
# For two-dimensional geometries
gmsh -2 file_path.geo -format msh2

# For three-dimensional geometries
gmsh -3 file_path.geo -format msh2
```

Convert `.MSH` file to `.XML`
```bash
dolfin-convert file_path.msh file_path.xml
```

To test if the resulting `.XML` file is a valid mesh for Dolfin try
```python
import dolfin
file_path = ... # Put the path to the .XML file here
mesh = dolfin.Mesh("file_path.xml")
dolfin.plot(mesh)
```

### 2D Case


### 3D Case


## Demo usage with FEniCS

Show how a domain in 2d can be created and a PDE can be solved on it using FEniCS.


Steps:
- Install Docker Desktop
- Pull the fenics docker image: ```docker pull quay.io/fenicsproject/stable:latest```
- Run the image:
    - docker run -ti -v $(pwd):/home/fenics/shared quay.io/fenicsproject/stable
    - docker run -ti quay.io/fenicsproject/stable:latest
