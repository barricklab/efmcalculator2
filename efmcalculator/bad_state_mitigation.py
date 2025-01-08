def detect_special_cases(sequence, circular=True):
    def can_tile(s):
        doubled_s = s + s
        modified_doubled = doubled_s[1:-1]
        return s in modified_doubled
    if can_tile(str(sequence).lower()) and circular:
        raise BadSequenceError("Circular sequence is an infinitely long SSR. Did you mean to use a linear strategy?")
    if len(str(sequence).lower()) < 21 and circular:
        raise BadSequenceError("EFM Calculator is limited to circular sequences of at least 21 bases.")


class BadSequenceError(ValueError):
    pass
