Install
====================================================

The gulls Simulator, a comprehensive tool built with C++, offers a robust platform for simulating microlensing events. This guide is tailored for users who wish to install the simulator without delving into the intricacies of its C++ codebase.

.. important::

   The gulls Simulator is designed for cross-platform compatibility. However, specific steps, such as installing dependencies, may vary slightly across different operating systems.


Installation Steps
------------------

1. **Essential Libraries**:

   Install essential C++ libraries used by gulls:

   - **GNU Scientific Library (GSL)**:
   
     .. code-block:: shell

        conda install -c conda-forge gsl

   - **CFITSIO** for FITS file handling:

     .. code-block:: shell

        conda install -c conda-forge cfitsio

2. **gulls Base Directory**:

   Define the base directory for gulls installation:

   .. code-block:: shell

      export GULLS_BASE_DIR='/path/to/gulls/'
      echo 'export GULLS_BASE_DIR="/path/to/gulls/"' >> ~/.bashrc  # For bash users

3. **Configuration File**:

   Although optional, creating a configuration file is advised:

   .. code-block:: shell

      touch ~/.gulls

4. **Compilation**:

   Navigate to the gulls root directory and compile the project:

   .. code-block:: shell

      cd $GULLS_BASE_DIR
      ./configure.sh
      make

   This process compiles the gulls Simulator, making it ready for use. Ensure your environment is correctly set up with a C++ compiler and necessary build tools.

.. note::

   The exact commands and steps, especially for environment and dependency setup, might need adjustments based on your operating system and the specific versions of software you're using. Always consult the official documentation of each software component for the most accurate guidance.

.. warning::

   The successful compilation of gulls depends on a correctly configured C++ development environment. If you encounter any issues during the `make` process, verify that all dependencies are correctly installed and that your C++ compiler supports the C++ standards required by gulls.
