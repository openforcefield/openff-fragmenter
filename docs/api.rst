API
===

Below is an outline of the API for the main functions of ``fragmenter`` See :doc:`examples <examples>` for details on
how to use these functions.

.. warning:: ``fragmenter`` The ``fragmenter`` package is still pre-alpha so the API is still in flux.
             pre-alpha so the API is still in flux.

Fragmentation Engines
---------------------

.. autopydantic_model:: fragmenter.fragment.Fragmenter
   :members:
   :undoc-members:

WBO
"""

.. autopydantic_model:: fragmenter.fragment.WBOFragmenter
   :members:
   :undoc-members:

.. autopydantic_model:: fragmenter.fragment.WBOOptions
   :members:
   :undoc-members:

Pfizer
""""""

.. autopydantic_model:: fragmenter.fragment.PfizerFragmenter
   :members:
   :undoc-members:

Fragmentation Outputs
---------------------

.. autopydantic_model:: fragmenter.fragment.Fragment
   :members:
   :undoc-members:

.. autopydantic_model:: fragmenter.fragment.FragmentationResult
   :members:
   :undoc-members:
