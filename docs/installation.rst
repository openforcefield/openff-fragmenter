Installing fragmenter
=====================

Installing using conda
----------------------
To install ``fragmenter`` with conda run the following

.. code-block:: bash

    conda install -c omnia fragmenter

Installing from source
----------------------
To install ``fragmenter`` from source, clone or download the `GitHub repo <https://github.com/openforcefield/fragmenter>`_.
From inside the ``fragmenter`` directory, run the following:

.. code-block:: bash

    python setup.py install

This command will not install dependencies. All dependencies are in the ``meta.yaml`` `file <https://github.com/openforcefield/fragmenter/blob/master/devtools/conda-envs/meta.yaml>`_

Prerequisites
-------------
``fragmenter`` is tested with python 3.6.

This toolkit uses `OpenEye <https://www.eyesopen.com/>`_ as a dependency so licenses for `oechem`, `oequacpac` and `oeomega` are required.

Warning
-------
``fragmenter`` is still pre-alpha. It is not fully tested and the API is still in flux.

