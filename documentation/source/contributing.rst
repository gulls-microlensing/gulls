Contributing to GULLS
=====================

Thank you for your interest in contributing to GULLS! This document provides guidelines for contributing to the project.

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

Building GULLS
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

- GULLS version/commit
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
