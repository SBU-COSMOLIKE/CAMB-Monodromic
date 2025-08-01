.. _fortran-compilers:

Fortran compilers
=================

CAMB internally uses modern (object-oriented) Fortran 2008 for most numerical calculations (see `docs <https://camb.info/doc/>`_),
and needs a fortran compiler to build the numerical library. The recommended compilers are

- gfortran
- Intel Fortran (ifort), version 18.0.1 or higher

The gfortran compiler is part of the standard "gcc" compiler package, and may be pre-installed on recent unix systems.
Check the version using "gfortran --version".

If you do not have a suitable Fortran compiler, you can get one as follows:

:Mac:
    Download the `binary installation <https://gcc.gnu.org/wiki/GFortranBinaries>`_
:Windows:
    Download gfortran as part of `MinGW-w64 <https://sourceforge.net/projects/mingw-w64/files>`_ (select x86_64 option in the installation program)
    or get latest from niXman on `GitHub <https://github.com/niXman/mingw-builds-binaries/releases>`_ (e.g. x86_64-13.2.0-release-win32-seh-msvcrt-rt_v11-rev1)
:Linux:
    To install from the standard repository use:

     - "sudo apt-get update; sudo apt-get install gfortran"

Alternatively you can compile and run in a container or virtual machine: e.g., see `CosmoBox <https://cosmologist.info/CosmoBox>`_.
For example, to run a configured shell in docker where you can install and run camb from the command line (after changing to the camb directory)::

    docker run -v /local/git/path/CAMB:/camb -i -t cmbant/cosmobox

Updating modified Fortran code
------------------------------

In the main CAMB source root directory, to re-build the Fortran binary including any
pulled or local changes use::

    python setup.py make

This will also work on Windows as long as you have MinGW-w64 installed under Program Files as described above.

NOTE: gfortran occasionally produces `memory leaks <https://gcc.gnu.org/bugzilla/show_bug.cgi?id=120637>`_: if you see
leaks running your code, try adding final methods to explicitly free any allocatable arrays or subcomponents.


Note that you will need to close all python instances using camb before you can re-load with an updated library.
This includes in Jupyter notebooks; just re-start the kernel or use::

    import IPython
    IPython.Application.instance().kernel.do_shutdown(True)

If you want to automatically rebuild the library from Jupyter you can do something like this::

    import subprocess
    import sys
    import os
    src_dir = '/path/to/git/CAMB'
    try:
        subprocess.check_output(r'python "%s" make'%os.path.join(src_dir, 'setup.py'),
                                stderr=subprocess.STDOUT)
        sys.path.insert(0,src_dir)
        import camb
        print('Using CAMB %s installed at %s'%(camb.__version__,
                                        os.path.dirname(camb.__file__)))

    except subprocess.CalledProcessError as E:
        print(E.output.decode())
