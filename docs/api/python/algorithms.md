---
og:description: "Every reconstruction algorithm in tigre.algorithms for Python, from FDK to SART, CGLS, LSQR and total-variation methods, with their options."
myst:
  html_meta:
    description: "Every reconstruction algorithm in tigre.algorithms for Python, from FDK to SART, CGLS, LSQR and total-variation methods, with their options."
---

# Python algorithms

In `tigre.algorithms`. All iterative algorithms are called as
`algorithm(proj, geo, angles, niter, **options)`; see
[](../../user-guide/algorithms.md) for which to choose.

```{eval-rst}
.. currentmodule:: tigre.algorithms

.. rubric:: Filtered backprojection

.. autosummary::
   :nosignatures:

   fdk
   fbp

.. rubric:: SART family

.. autosummary::
   :nosignatures:

   sart
   sirt
   ossart
   sart_tv
   ossart_tv

.. rubric:: Krylov subspace methods

.. autosummary::
   :nosignatures:

   cgls
   lsqr
   lsmr
   hybrid_lsqr
   ab_gmres
   ba_gmres
   irn_tv_cgls
   irn_piccs_cgls
   hybrid_flsqr_tv

.. rubric:: Total variation (POCS family)

.. autosummary::
   :nosignatures:

   asd_pocs
   b_asd_pocs_beta
   os_asd_pocs
   awasd_pocs
   os_awasd_pocs
   pcsd
   os_pcsd
   aw_pcsd
   os_aw_pcsd

.. rubric:: Prior image (PICCS)

.. autosummary::
   :nosignatures:

   piccs
   os_piccs

.. rubric:: Proximal methods

.. autosummary::
   :nosignatures:

   fista
   ista

.. rubric:: Statistical methods

.. autosummary::
   :nosignatures:

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

.. autofunction:: irn_piccs_cgls

.. autofunction:: hybrid_flsqr_tv

Total variation (POCS family)
-----------------------------

.. autofunction:: asd_pocs

.. autofunction:: b_asd_pocs_beta

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
