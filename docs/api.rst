API
===

Below is an outline of the API for the main functions of ``openff-fragmenter`` See the examples for details on how to
use these objects.

.. warning:: The ``openff-fragmenter`` package is still pre-alpha so the API is still in flux.

Fragmentation Engines
---------------------

.. autopydantic_model:: openff.fragmenter.fragment.Fragmenter
   :members:
   :undoc-members:

WBO
"""

.. autopydantic_model:: openff.fragmenter.fragment.WBOFragmenter
   :members:
   :undoc-members:

.. autopydantic_model:: openff.fragmenter.fragment.WBOOptions
   :members:
   :undoc-members:

Pfizer
""""""

.. autopydantic_model:: openff.fragmenter.fragment.PfizerFragmenter
   :members:
   :undoc-members:

Fragmentation Outputs
---------------------

.. autopydantic_model:: openff.fragmenter.fragment.Fragment
   :members:
   :undoc-members:

.. autopydantic_model:: openff.fragmenter.fragment.FragmentationResult
   :members:
   :undoc-members:
