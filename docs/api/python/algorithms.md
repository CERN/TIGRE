# Algorithms

In `tigre.algorithms`. All iterative algorithms are called as
`algorithm(proj, geo, angles, niter, **options)`; see
[](../../user-guide/algorithms.md) for which to choose.

```{eval-rst}
.. currentmodule:: tigre.algorithms

.. autosummary::
   :nosignatures:

   fdk
   fbp
   sart
   sirt
   ossart
   sart_tv
   ossart_tv
   cgls
   lsqr
   lsmr
   hybrid_lsqr
   ab_gmres
   ba_gmres
   irn_tv_cgls
   hybrid_flsqr_tv
   asd_pocs
   os_asd_pocs
   awasd_pocs
   os_awasd_pocs
   pcsd
   os_pcsd
   aw_pcsd
   os_aw_pcsd
   piccs
   os_piccs
   fista
   ista
   mlem

Filtered backprojection
-----------------------

.. autofunction:: fdk

.. autofunction:: fbp

SART family
-----------

.. autofunction:: sart

.. autofunction:: sirt

.. autofunction:: ossart

.. autofunction:: sart_tv

.. autofunction:: ossart_tv

Krylov subspace methods
-----------------------

.. autofunction:: cgls

.. autofunction:: lsqr

.. autofunction:: lsmr

.. autofunction:: hybrid_lsqr

.. autofunction:: ab_gmres

.. autofunction:: ba_gmres

.. autofunction:: irn_tv_cgls

.. autofunction:: hybrid_flsqr_tv

Total variation (POCS family)
-----------------------------

.. autofunction:: asd_pocs

.. autofunction:: os_asd_pocs

.. autofunction:: awasd_pocs

.. autofunction:: os_awasd_pocs

.. autofunction:: pcsd

.. autofunction:: os_pcsd

.. autofunction:: aw_pcsd

.. autofunction:: os_aw_pcsd

Prior image (PICCS)
-------------------

.. autofunction:: piccs

.. autofunction:: os_piccs

Proximal methods
----------------

.. autofunction:: fista

.. autofunction:: ista

Statistical methods
-------------------

.. autofunction:: mlem

```
