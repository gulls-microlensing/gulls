Quickstart
============
.. raw:: html

   <hr>

.. _general_quickstart:

Getting Started
---------------

1. **Installation**: Follow the instructions in :doc:`Installation <install_gulls>`
2. **Prepare Input Files**: Create your parameter file and catalogs (see :doc:`Input Formats <input_formats>`)
3. **Validate Inputs** (Recommended):

   .. code-block:: bash
   
      python3 scripts/validate_inputs.py your_parameter_file.prm

   This checks for common catalog issues before running simulations.

4. **Run Simulation**: See :doc:`Running Gulls <run_overview>` for execution details

.. tip::
   Start with the smoke test examples in ``smoke_test/parameterfiles/`` to understand
   the expected input format.
