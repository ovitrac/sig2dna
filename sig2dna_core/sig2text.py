"""
sig2text — typed, evidential serialization of chemical measurements
(EXPERIMENTAL: schema under active development; the field set may change
between 0.x versions).

The design goal: compile deterministic measurements into a *chemical
document* an LLM can reason over safely. The compiler is deterministic;
the LLM only ever renders prose FROM documents (Level 2) and never
produces the measurements (Level 0) or the typed claims (Level 1).

Principles encoded here:

- **Typed namespaces** so lexical collisions are impossible: ``sigma:``
  signal morphology, ``eps:`` events/channels, ``phi:`` phenotype state,
  ``mu:`` mechanism/interpretation.
- **Evidentiality as a type/scope system**: every claim carries state,
  magnitude, validity scope, and evidence status (measured / attributed /
  candidate / suppressed). One phenomenon may have several
  representations (level, profile, sign) with DIFFERENT transport
  scopes; the language never forces them to share one.
- **Masks**: missing is never zero; a masked claim keeps its reason.
- **STATE and CAUSE never mix**: measured states live in ``phi:``/
  ``eps:``; interpretations live in ``mu:`` with an explicit basis.
- **Suppressed negative findings** are serialized as such -- a validated
  "not identifiable" is information, not absence.
- **Round-trip law**: ``parse(serialize(doc)) == doc`` for the scientific
  state object; serialization must not lose or coerce.
- **Calibration is part of the claim**: a score carries the mode that
  produced it (``replicate`` / ``reference_population`` /
  ``density_conditioned`` / ``global_proxy``, in decreasing evidential
  quality); a density-conditioned score must never masquerade as
  replicate-grade significance.
"""
from __future__ import annotations

from dataclasses import dataclass, field

SCHEMA = "0.1-experimental"
NAMESPACES = ("sigma", "eps", "phi", "mu")
STATUSES = ("measured", "attributed", "candidate", "suppressed")


@dataclass(frozen=True)
class Claim:
    """One typed claim. ``value=None`` with ``mask_reason`` = masked
    (missing != zero); ``status='suppressed'`` = validated negative
    finding whose token is deliberately withheld."""
    namespace: str
    name: str
    state: str = ""
    value: float | None = None
    scope: str = ""
    status: str = "measured"
    mask_reason: str = ""
    calibration: str = ""

    def __post_init__(self):
        if self.namespace not in NAMESPACES:
            raise ValueError(f"unknown namespace {self.namespace!r}")
        if self.status not in STATUSES:
            raise ValueError(f"unknown status {self.status!r}")
        if self.value is None and self.status == "measured" \
                and not self.mask_reason:
            raise ValueError("missing value requires a mask_reason "
                             "(missing is never zero)")


@dataclass
class Document:
    sample_id: str
    domain: str
    claims: list = field(default_factory=list)

    def serialize(self) -> str:
        lines = [f"SIG2TEXT schema={SCHEMA}",
                 f"  sample_id: {self.sample_id}",
                 f"  domain: {self.domain}"]
        for c in self.claims:
            if c.value is None:
                body = f"value=NA mask_reason={c.mask_reason!r}"
            else:
                body = f"value={c.value!r}"
            cal = f" calibration={c.calibration}" if c.calibration else ""
            lines.append(f"  {c.namespace}:{c.name} state={c.state!r} "
                         f"{body} scope={c.scope!r} status={c.status}{cal}")
        return "\n".join(lines) + "\n"


def parse(text: str) -> Document:
    """Inverse of :meth:`Document.serialize` (the round-trip law)."""
    import re
    lines = text.strip().split("\n")
    if not lines[0].startswith("SIG2TEXT"):
        raise ValueError("not a sig2text document")
    doc = Document(sample_id=lines[1].split(": ", 1)[1],
                   domain=lines[2].split(": ", 1)[1])
    pat = re.compile(
        r"^\s*(\w+):(\S+) state='([^']*)' "
        r"(?:value=(\S+)|value=NA mask_reason='([^']*)') "
        r"scope='([^']*)' status=(\w+)(?: calibration=(\S+))?$")
    for line in lines[3:]:
        m = pat.match(line)
        if not m:
            raise ValueError(f"unparseable claim line: {line!r}")
        ns, name, state, val, mask, scope, status, cal = m.groups()
        doc.claims.append(Claim(
            namespace=ns, name=name, state=state,
            value=None if val is None else float(val),
            scope=scope, status=status,
            mask_reason=mask or "", calibration=cal or ""))
    return doc
