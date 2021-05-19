=================
OpenFF Fragmenter
=================

The main purpose of ``openff-fragmenter`` is to fragment molecules for quantum chemical (QC) torsion drives.

.. warning:: ``openff-fragmenter`` is still pre-alpha. It is not fully tested and the API is still in flux.

Currently two fragmentation schemes are supported:

.. rst-class:: spaced-list

 * a Wiberg Bond Order (WBO) sensitive fragmentation scheme :cite:`2020:stern` *(recommended)*
 * the fragmentation schema detailed by Rai *et al* :cite:`2019:bakken` reffered to in this
   package as the 'Pfizer' scheme.

WBO Sensitive Fragmentation
---------------------------

The assumption when fragmenting molecules is that the chemistry is localized and that removing or changing remote
substituents (defined as substituents more than 2 bonds away from the central bond that is being driven in the torsion
drive) will not change the torsion potential around the bond of interest. However, that is not always the case.
``openff-fragmenter`` uses the Wiberg Bond Order (WBO) as a surrogate signal to determine if the chemistry around the bond of
interest was destroyed during fragmentation relative to the bond in the parent molecule.

The WBO is a measure of electronic population overlap between two atoms in a bond. It can be quickly calculated from
an empirical QC calculation and is given by:

.. math:: W_{AB} = \sum_{\mu\in{A}}\sum_{\nu\in{B}}|D_{\mu\nu}|^2

Where :math:`A` and :math:`B` are atoms :math:`A` and :math:`B` in a bond, :math:`D` is the density matrix and :math:`\mu`
and :math:`\nu` are occupied orbitals on atoms :math:`A` and :math:`B` respectively.

``openff-fragmenter`` calculates the WBO of the parent molecules, then fragments according to a set of rules and then
recalculates the WBO of the fragments. If the WBO for the fragment of the bond of interest changes more than a user's
specified threshold, ``openff-fragmenter`` will add more substituents until the WBO of the bond of interest is within the user
specified threshold.

Contributors
------------
* `Chaya D. Stern (MSKCC / Weill Cornell) <https://github.com/ChayaSt>`_
* `John D. Chodera (MSKCC) <https://github.com/jchodera>`_

Acknowledgment
--------------
CDS is funded by a fellowship from `The Molecular Sciences Software Institute <http://molssi.org/>`_

References
----------

.. bibliography:: index.bib
    :cited:
    :style: unsrt

.. toctree::
   :maxdepth: 2
   :hidden:

   Overview <self>
   installation
   Examples <fragment-molecules.ipynb>
   api