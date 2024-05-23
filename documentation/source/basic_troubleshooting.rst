Troubleshooting
======================================


For any issues or questions related to the GULLS Simulator, we strongly encourage the use of the `GitHub issue tracking system <https://github.com/gulls-microlensing/gulls/issues>`_. This platform streamlines the process of managing and resolving issues by keeping a clear record of active problems and their status.

.. figure:: ../build/html/_static/github_issueTracking.png
   :align: center
   :width: 700
   :alt: GitHub Issues Page

.. raw:: html

   <hr>

To submit a new issue, click the green ``New Issue`` button. This action will navigate you to a submission form, similar to what's depicted below:

.. figure:: ../build/html/_static/github_newIssue.png
   :align: center
   :width: 700
   :alt: GitHub New Issue Page

.. raw:: html

   <hr>

For an efficient and swift troubleshooting process, please provide the following information in your issue report:

* **Problem Description**: Offer a concise explanation of the issue.
* **Attempted Solutions**: Describe what you have tried in order to resolve the issue.
* **Required Files**: Attach the ``parameterfile.prm`` if applicable.
* **System Information** (for installation or execution issues):
  * Operating System
  * GCC version (obtain via ``gcc --version``) or Clang version (``clang --version``)

.. raw:: html

   <hr>

Common Issues and Solutions
---------------------------

**White Square Images**:
If your simulation output is a white square, this could indicate an issue with magnitude values in your star, source, or lens catalogs, or you might have encountered a super-bright star in the image. To diagnose:

* Ensure ``zscale`` is selected in DS9.
* Generate a larger ``PRETTY_PIC`` and inspect for super-bright stars, which appear as white squares. If prevalent, reviewing the color magnitude diagrams of your star lists is advisable.

**Debug Information**:
For additional insights, you can increase the level of debug information by repeating the ``-d`` flags when running the simulator.

**Valid Stars Error**:
Receiving a message about no valid stars suggests a misconfiguration of ``NFILTERS``. Verify this setting corresponds correctly to your simulation requirements.
