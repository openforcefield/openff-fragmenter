Installation
============

Installing using mamba
----------------------

The recommended way to install ``openff-fragmenter`` is via the ``mamba`` package manger:

.. code-block:: bash

    mamba install -c conda-forge openff-fragmenter

If you do not have Conda installed, see the `OpenFF installation guide <openff.docs:install>`_.

If you have access to the OpenEye toolkits (namely ``oechem``, ``oequacpac`` and ``oeomega``) we recommend installing
these also as these can speed up fragmentation times significantly:

.. code-block:: bash

    mamba install -c openeye openeye-toolkits

Installing from source
----------------------

To install ``openff-fragmenter`` from source begin by cloning the repository from `github
<https://github.com/openforcefield/openff-fragmenter>`_:

.. code-block:: bash

    git clone https://github.com/openforcefield/openff-fragmenter.git
    cd fragmenter

Create a custom conda environment which contains the required dependencies and activate it:

.. code-block:: bash

    mamba env create --name fragmenter --file devtools/conda-envs/meta.yaml
    mamba activate fragmenter

Finally, install ``openff-fragmenter`` itself:

.. code-block:: bash

    python setup.py develop
