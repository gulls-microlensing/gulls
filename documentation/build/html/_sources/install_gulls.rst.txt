Install
==========

The Gulls simulator, a comprehensive tool built with C++, offers a robust platform for simulating microlensing events. This guide is tailored for users who wish to install the simulator without delving into the intricacies of its C++ codebase.

.. important::

   The Gulls simulator is designed for cross-platform compatibility. However, specific steps, such as installing dependencies, may vary slightly across different operating systems.


**Prerequisites**

Before installing Gulls, ensure you have:

- **C++17** compatible compiler (GCC, Clang, or MSVC)
- **Fortran** compiler (gfortran recommended)
- **Python 3.7+** (for validation tools and smoke tests)

For the **recommended CMake build**:

- **CMake** 3.10 or later

For the **legacy Makefile build** (OS/environment dependent):

- Traditional build tools (make, etc.)

**Installation Methods**

Gulls supports two build systems. **We recommend using CMake** for most users as it handles
dependencies more reliably across different platforms.

Method 1: CMake (Recommended)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

1. **Install Dependencies via Conda** (recommended):

   The easiest way to install all dependencies is using the provided conda environment:

   .. code-block:: shell

      conda env create -f environment.yml
      conda activate gulls

   This installs:
   
   - GNU Scientific Library (GSL)
   - CFITSIO (for FITS file handling)
   - CMake and compilers
   - Python packages for validation and testing

   **Or install manually**:

   .. code-block:: shell

      conda install -c conda-forge gsl cfitsio cmake compilers

2. **Build with CMake**:

   Navigate to the Gulls root directory and compile:

   .. code-block:: shell

      cd /path/to/gulls/
      cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
      cmake --build build --parallel

   The executables will be in the ``bin/`` directory:
   
   - ``bin/gulls_std.x`` - Standard microlensing simulation
   - ``bin/gulls_croin.x`` - Crowding analysis
   - ``bin/gullsFish.x`` - Fisher matrix analysis

   .. tip::
      For development or debugging, use ``-DCMAKE_BUILD_TYPE=Debug`` instead of ``Release``.
      The ``--parallel`` flag speeds up compilation by using multiple CPU cores.

3. **Verify Installation**:

   Test your installation by running the smoke tests:

   .. code-block:: shell

      python3 smoke_test/run_smoke_test.py --ci

   This runs a subset of quick validation tests to ensure the build is working correctly.

Method 2: Legacy Makefile (Alternative)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. warning::
   This build method is OS and environment dependent and may require manual configuration.
   Use the CMake method above unless you have specific reasons to use the legacy build system.

1. **Install Dependencies**:

   Install libraries manually (commands vary by OS):

   .. code-block:: shell

      # Using conda
      conda install -c conda-forge gsl cfitsio

   .. code-block:: shell

      # Or using apt (Ubuntu/Debian)
      sudo apt-get install libgsl-dev libcfitsio-dev

   .. code-block:: shell

      # Or using homebrew (macOS)
      brew install gsl cfitsio

2. **Set Environment Variables**:

   Define the base directory for gulls installation:

   .. code-block:: shell

      export GULLS_BASE_DIR='/path/to/gulls/'
      echo 'export GULLS_BASE_DIR="/path/to/gulls/"' >> ~/.bashrc  # For bash users

3. **Configuration File** (optional):

   .. code-block:: shell

      touch ~/.gulls

4. **Compilation**:

   Navigate to the Gulls root directory and compile:

   .. code-block:: shell

      cd $GULLS_BASE_DIR
      ./configure.sh
      make

   The executables will be in the ``bin/`` directory.

.. note::

   The exact commands and steps, especially for environment and dependency setup, might need adjustments based on your operating system and the specific versions of software you're using. Always consult the official documentation of each software component for the most accurate guidance.

.. warning::

   The successful compilation of Gulls depends on a correctly configured C++ development environment. If you encounter any issues during the `make` process, verify that all dependencies are correctly installed and that your C++ compiler supports the C++ standards required by Gulls.
