API
===

Below is an outline of the API for the main functions of ``fragmenter`` See :doc:`examples <examples>` for details on
how to use these functions.

.. warning:: ``fragmenter`` The ``fragmenter`` package is still pre-alpha so the API is still in flux.is still pre-alpha so the API is still in flux.

Fragmentation Engines
---------------------

.. currentmodule:: fragmenter.fragment

.. autoclass:: WBOOptions
    :members:

.. autoclass:: WBOFragmenter
    :members: +fragment, json, dict

.. autoclass:: PfizerFragmenter
    :members: +fragment, json, dict

.. autoclass:: Fragment
    :members:

.. autoclass:: FragmentationResult
    :members:
