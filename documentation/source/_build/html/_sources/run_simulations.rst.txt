Running gulls
===============================

After configuring your environment and meeting all prerequisites, follow these steps to run simulations with the GULLS Simulator.

Building the Simulator
----------------------

1. **Open the Makefile**:
   Navigate to the GULLS source directory and open the ``Makefile``:

   .. code-block:: shell

      cd /path/to/gulls/src/
      open Makefile

2. **Configure the Makefile**:
   Edit the Makefile to set the ``BASEDIR`` and compiler flags such as ``FFLAGS``, ``CFLAGS``, and ``LINKERFLAGS`` to match your system's configuration.

3. **Select an Executable**:
   Based on your simulation needs, choose an executable from the list defined in the Makefile. For instance, to simulate single lens microlensing events, you would focus on ``gullsSingle``.

4. **Build the Executable**:
   Compile the chosen executable by running:

   .. code-block:: shell

      make gullsSingle

   Clean up any intermediate files if necessary:

   .. code-block:: shell

      make clean


Running a Simulation
--------------------

1. **Change to the Bin Directory**:
   Once the build process completes, navigate to the bin directory:

   .. code-block:: shell

      cd /path/to/gulls/bin

2. **Execute the Simulator**:
   Launch the GULLS simulator by running the built executable. For example:

   .. code-block:: shell

      ./gullsSingle.x

   You'll be prompted with usage guidelines. Follow these to run your simulation effectively.

3. **Simulation Input**:
   Prepare your input parameter file (e.g., ``singlelens.prm``) and specify it when running the simulator. Adjust the flags according to your needs:

   - **-i**: Path to the input parameter file.
   - **-s**: Simulation instance number, usually set to 0 for initial runs.
   - **-f**: Field number for the observation, if applicable.
   - **-d**: Debug mode, can be repeated for increased verbosity.

   Example command:

   .. code-block:: shell

      ./gullsSingle.x -i PATH/singlelens.prm -s 0 -f 83 -d

   This command initiates the simulation with specified parameters and outputs data to the ``OUTPUT_DIR`` configured in your parameter file.

Congratulations! You've successfully run a simulation with the GULLS Simulator. Explore the output in your specified directory, and adjust parameters for further simulations as needed.
