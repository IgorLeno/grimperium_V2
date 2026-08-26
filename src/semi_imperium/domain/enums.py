"""Explicit lifecycle vocabulary for Semi-Imperium calculations.

Every state a molecule can be in is a named enum member. There is no
"empty means pending" or "null means not verified" convention here: a
record that has not been verified says ``NOT_REQUESTED`` or ``PENDING``
explicitly, so a reader never has to guess what a blank cell meant.
"""

from __future__ import annotations

from enum import Enum


class CalculationState(str, Enum):
    """Lifecycle state of one individual calculation.

    ``VERIFIED`` and ``UNVERIFIED`` both mean the calculation produced a
    result; they differ only in whether the minimum-verification policy
    confirmed the geometry is a true minimum. ``SADDLE`` is a scientific
    outcome, not an error: the calculation succeeded and the verification
    proved the geometry is a saddle point.
    """

    PENDING = "pending"
    """Registered and waiting to run. Nothing has been executed yet."""

    RUNNING = "running"
    """Execution started and has not reached a terminal outcome."""

    VERIFIED = "verified"
    """Finished with a result whose geometry was confirmed as a minimum."""

    UNVERIFIED = "unverified"
    """Finished with a result, but minimum verification was not confirmed."""

    SADDLE = "saddle"
    """Finished with a result whose geometry was proven to be a saddle."""

    FAILED = "failed"
    """Execution ended without a usable scientific result."""

    @property
    def is_terminal(self) -> bool:
        """Whether no further transition is allowed from this state."""
        return self in TERMINAL_CALCULATION_STATES

    @property
    def has_result(self) -> bool:
        """Whether the state carries a scientific result worth persisting."""
        return self in RESULT_CALCULATION_STATES


class VerificationPolicy(str, Enum):
    """How strongly a run demands proof that a geometry is a minimum.

    The policy is part of the calculation signature, so tightening it
    invalidates reuse of results produced under a looser policy.
    """

    NONE = "none"
    """No frequency verification is performed."""

    ADVISORY = "advisory"
    """Frequencies are computed and recorded, but never reject a result."""

    REQUIRE_MINIMUM = "require_minimum"
    """A result is only accepted when all frequencies are real."""


class VerificationOutcome(str, Enum):
    """Result of applying the verification policy to one calculation."""

    NOT_REQUESTED = "not_requested"
    """Policy is ``NONE``; verification was deliberately not attempted."""

    PENDING = "pending"
    """Policy requires verification and it has not produced a verdict yet."""

    MINIMUM_CONFIRMED = "minimum_confirmed"
    """All frequencies are real: the geometry is a minimum."""

    SADDLE_POINT = "saddle_point"
    """At least one imaginary frequency: the geometry is a saddle point."""

    INCONCLUSIVE = "inconclusive"
    """Verification ran but could not decide (e.g. frequencies unreadable)."""

    FAILED = "failed"
    """Verification itself failed to execute."""


class RunState(str, Enum):
    """Lifecycle state of a run, i.e. one batch of calculations."""

    PENDING = "pending"
    RUNNING = "running"
    COMPLETED = "completed"
    PARTIAL = "partial"
    FAILED = "failed"
    CANCELLED = "cancelled"

    @property
    def is_terminal(self) -> bool:
        """Whether no further transition is allowed from this state."""
        return self in TERMINAL_RUN_STATES


TERMINAL_CALCULATION_STATES = frozenset(
    {
        CalculationState.VERIFIED,
        CalculationState.UNVERIFIED,
        CalculationState.SADDLE,
        CalculationState.FAILED,
    }
)

RESULT_CALCULATION_STATES = frozenset(
    {
        CalculationState.VERIFIED,
        CalculationState.UNVERIFIED,
        CalculationState.SADDLE,
    }
)

TERMINAL_RUN_STATES = frozenset(
    {
        RunState.COMPLETED,
        RunState.PARTIAL,
        RunState.FAILED,
        RunState.CANCELLED,
    }
)

#: Explicit transition matrix. Terminal states accept no further transitions.
ALLOWED_CALCULATION_TRANSITIONS: dict[CalculationState, frozenset[CalculationState]] = {
    CalculationState.PENDING: frozenset(
        {CalculationState.RUNNING, CalculationState.FAILED}
    ),
    CalculationState.RUNNING: frozenset(
        {
            CalculationState.VERIFIED,
            CalculationState.UNVERIFIED,
            CalculationState.SADDLE,
            CalculationState.FAILED,
        }
    ),
    CalculationState.VERIFIED: frozenset(),
    CalculationState.UNVERIFIED: frozenset(),
    CalculationState.SADDLE: frozenset(),
    CalculationState.FAILED: frozenset(),
}

ALLOWED_RUN_TRANSITIONS: dict[RunState, frozenset[RunState]] = {
    RunState.PENDING: frozenset(
        {RunState.RUNNING, RunState.CANCELLED, RunState.FAILED}
    ),
    RunState.RUNNING: frozenset(
        {
            RunState.COMPLETED,
            RunState.PARTIAL,
            RunState.FAILED,
            RunState.CANCELLED,
        }
    ),
    RunState.COMPLETED: frozenset(),
    RunState.PARTIAL: frozenset(),
    RunState.FAILED: frozenset(),
    RunState.CANCELLED: frozenset(),
}

#: Verification outcomes that are coherent with each calculation state.
#:
#: This is what keeps "verified" honest: a record cannot claim
#: ``VERIFIED`` while carrying a ``SADDLE_POINT`` verdict.
COHERENT_VERIFICATION_OUTCOMES: dict[
    CalculationState, frozenset[VerificationOutcome]
] = {
    CalculationState.PENDING: frozenset(
        {VerificationOutcome.NOT_REQUESTED, VerificationOutcome.PENDING}
    ),
    CalculationState.RUNNING: frozenset(
        {VerificationOutcome.NOT_REQUESTED, VerificationOutcome.PENDING}
    ),
    CalculationState.VERIFIED: frozenset({VerificationOutcome.MINIMUM_CONFIRMED}),
    CalculationState.UNVERIFIED: frozenset(
        {
            VerificationOutcome.NOT_REQUESTED,
            VerificationOutcome.INCONCLUSIVE,
            VerificationOutcome.FAILED,
        }
    ),
    CalculationState.SADDLE: frozenset({VerificationOutcome.SADDLE_POINT}),
    CalculationState.FAILED: frozenset(
        {
            VerificationOutcome.NOT_REQUESTED,
            VerificationOutcome.PENDING,
            VerificationOutcome.FAILED,
            VerificationOutcome.INCONCLUSIVE,
        }
    ),
}

#: States a stored calculation must be in before it can be reused.
#:
#: ``SADDLE`` is included because it is a reproducible scientific finding:
#: rerunning the same signature on the same molecule would find the same
#: saddle point. Callers that need a true minimum pass a narrower set to
#: :meth:`semi_imperium.persistence.SemiImperiumStore.find_reusable`.
DEFAULT_REUSABLE_STATES = frozenset(RESULT_CALCULATION_STATES)

#: Narrower reuse set for callers that only accept confirmed minima.
VERIFIED_ONLY_REUSABLE_STATES = frozenset({CalculationState.VERIFIED})


__all__ = [
    "ALLOWED_CALCULATION_TRANSITIONS",
    "ALLOWED_RUN_TRANSITIONS",
    "COHERENT_VERIFICATION_OUTCOMES",
    "DEFAULT_REUSABLE_STATES",
    "RESULT_CALCULATION_STATES",
    "TERMINAL_CALCULATION_STATES",
    "TERMINAL_RUN_STATES",
    "VERIFIED_ONLY_REUSABLE_STATES",
    "CalculationState",
    "RunState",
    "VerificationOutcome",
    "VerificationPolicy",
]
