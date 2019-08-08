Fragmenter API
==============

``Fragmenter`` is pre-alpha so the API is still in flux.

Below is an outline of the API for the main functions of ``fragmenter``
See :doc:`examples <../examples>` for details on how to use these functions.

enumerate_states
----------------

Ionization can have a profound impact on the electronic structure and the torsion barrier which can lead to different
fragmentation decisions. Therefore, it is important to enumerate reasonable ionization states at physiological pH. ``fragmenter``
provides a wrapper around OpenEye's

.. toctree::
    :maxdepth: 2

    states.rst


fragment
--------

The `fragment` module is the core of ``fragmenter``. It provides two ways to fragment molecules:

1. Combinatorial fragmentation
    This scheme generates all possible fragments for a molecule without fragmenting rings and selected functional groups.
    It is not recommended for general use. It was used to generate the benchmark set used to validate ``fragmenter``
2. WBO fragmentation
    This scheme uses the change in WBO in the rotatable bonds to decide if the fragment needs to continue being grown out.
    The threshold for this change can be provided by the user. The default and recommended threshold is 0.01.

.. toctree::
    :maxdepth: 2

   fragment.rst





