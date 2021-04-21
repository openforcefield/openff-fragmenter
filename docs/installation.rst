Installation
============

Installing using conda
----------------------

The recommended way to install ``fragmenter`` is via the ``conda`` package manger:

.. code-block:: bash

    conda install -c omnia fragmenter

If you have access to the OpenEye toolkits (namely ``oechem``, ``oequacpac`` and ``oeomega``) we recommend installing
these also as these can speed up fragmentation times significantly:

.. code-block:: bash

    conda install -c openeye openeye-toolkits

Installing from source
----------------------

To install ``fragmenter`` from source begin by cloning the repository from `github
<https://github.com/openforcefield/fragmenter>`_:

.. code-block:: bash

    git clone https://github.com/openforcefield/fragmenter.git
    cd fragmenter

Create a custom conda environment which contains the required dependencies and activate it:

.. code-block:: bash

    conda env create --name fragmenter --file devtools/conda-envs/meta.yaml
    conda activate fragmenter

Finally, install ``fragmenter`` itself:

.. code-block:: bash

    python setup.py develop
