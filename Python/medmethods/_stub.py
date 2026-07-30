"""Helper for the not-yet-ported methods."""
from __future__ import annotations


def make_stub(name, title, r_fn, reason):
    """Build a function that raises NotImplementedError with a useful pointer."""
    def _stub(*args, **kwargs):
        raise NotImplementedError(
            "%s (%s) is not yet ported to Python. %s "
            "Use the R implementation MedMethods::%s() meanwhile; see the "
            "package README for the porting roadmap." % (name, title, reason, r_fn)
        )
    _stub.__name__ = name
    _stub.__doc__ = "%s -- not yet ported. %s Use R's MedMethods::%s()." % (
        title, reason, r_fn)
    return _stub
