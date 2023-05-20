"""Test mesh conversion.

Gmsh meshes are saved as .MSH files, these can be converted to dolfin-readable
.XML files using the `dolfin-convert` command. Depending on the Gmsh version the
conversion might cause problems.

With Gmsh version 4.11.1 create the .MSH file with the `-format msh2` option, i.e.
`gmsh -2 file.geo -format msh2` or `gmsh -3 file.geo -format msh2`
"""
if __name__ == "__main__":
    import os
    from pathlib import Path

    import dolfin

    ROOT_DIR = Path("/home/fenics/shared/")
    TEST_FILES = [
        ROOT_DIR / "tests/data/2d_domain_with_hole.msh",
        ROOT_DIR / "tests/data/sphere_3d.msh",
    ]

    # Convert .MSH to .XML
    for fn in TEST_FILES:
        xml_file = fn.with_suffix(".xml")
        os.system(f"dolfin-convert {fn} {xml_file}")

    # Read mesh into dolfin and plot it
    for fn in TEST_FILES:
        dolfin.Mesh(str(fn.with_suffix(".xml")))
