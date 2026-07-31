"""Constraint vocabulary for the constraints package.

The :class:`~bt4.domain.contracts.Constraint` protocol itself lives in ``domain``
(so the optimizer can consume it while every pure layer imports only ``domain``);
it is re-exported here so ``from bt4.constraints import Constraint`` keeps
working.
"""

from __future__ import annotations

from bt4.domain.contracts import Constraint

__all__ = ["Constraint"]
