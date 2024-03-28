Troubleshooting
======================================


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
