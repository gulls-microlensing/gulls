Contributing to Gulls
=====================

Thank you for your interest in contributing to Gulls! This document provides guidelines for contributing to the project.

Getting Started
---------------

1. Fork the repository
2. Clone your fork: ``git clone <your-fork-url>``
3. Create a feature branch: ``git checkout -b feature/your-feature-name``
4. Make your changes
5. Test your changes: ``python3 smoke_test/run_smoke_test.py``
6. Submit a pull request

Development Workflow
--------------------

Building Gulls
~~~~~~~~~~~~~~

.. code-block:: bash

   cmake -S . -B build -DCMAKE_BUILD_TYPE=Debug
   cmake --build build

Running Tests
~~~~~~~~~~~~~

.. code-block:: bash

   # Run all smoke tests
   python3 smoke_test/run_smoke_test.py

   # Run CI subset (faster)
   python3 smoke_test/run_smoke_test.py --ci

   # Run specific test
   python3 smoke_test/run_smoke_test.py --cases std-binary

Code Style
~~~~~~~~~~

- Follow existing C++ style conventions
- Use meaningful variable and function names
- Add comments for complex algorithms
- Update documentation for new features

Types of Contributions
----------------------

Bug Reports
~~~~~~~~~~~

When reporting bugs, please include:

- Gulls version/commit
- Operating system
- Steps to reproduce
- Expected vs actual behavior
- Relevant log files

Feature Requests
~~~~~~~~~~~~~~~~

For new features, please:

- Check existing issues first
- Provide a clear description
- Explain the use case
- Consider backward compatibility

Code Contributions
~~~~~~~~~~~~~~~~~~

- Keep changes focused and atomic
- Add tests for new functionality
- Update documentation
- Ensure all tests pass

Testing
-------

All contributions must pass the smoke tests:

.. code-block:: bash

   python3 smoke_test/run_smoke_test.py

The CI system will automatically run tests on pull requests.

Adding Validation for New Error Conditions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If you encounter a condition that causes a **breaking computational or science error**,
please add validation to catch it early:

**When to Add Validation:**

- ✅ Conditions that cause crashes, hangs, or incorrect physics
- ✅ Invalid catalog configurations (e.g., all sources closer than all lenses)
- ✅ Missing required columns for specific simulation modes
- ✅ Numerical values that cause undefined behavior (NaN, Inf, negative distances)
- ❌ Warnings or soft errors that don't break the simulation
- ❌ Performance issues that don't affect correctness

**How to Add Validation:**

1. Add a validation function to ``smoke_test/validation.py``:

   .. code-block:: python

      def verify_your_condition(params: Dict[str, str]) -> None:
          """Check that [condition] is satisfied."""
          # Load and check relevant data
          if condition_violated:
              raise SmokeTestError(
                  "Clear error message explaining what's wrong "
                  "and how to fix it"
              )

2. Export it in ``__all__`` at the bottom of ``validation.py``

3. Call it in ``smoke_test/runner.py`` in the validation loop (around line 141)

4. The validation will automatically be used by:
   
   - CI smoke tests
   - ``scripts/validate_inputs.py`` (for users)

**Example:** The infinite loop bug we encountered could have been caught by validating
that at least some source/lens distance pairs are valid (source > lens). This is now
implemented in ``verify_source_lens_compatibility()``.

.. tip::
   Write clear, actionable error messages. Users should understand what's wrong
   and how to fix it without reading the code.

Documentation
-------------

When adding new features:

- Update relevant documentation files
- Add examples if appropriate
- Update parameter reference if new parameters are added

Pull Request Process
--------------------

1. Ensure your branch is up to date with main
2. Run smoke tests locally
3. Create a clear, descriptive pull request
4. Reference any related issues
5. Respond to review feedback promptly

Questions?
----------

Feel free to open an issue for questions about contributing or the codebase.
