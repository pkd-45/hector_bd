"""Regression for scripts/prepare_pilot_inputs.py sample comparison (Claude review of 489e737, 11 Nov 2026).
R's jsonlite serialises doubles with 15 significant digits; the comparison must accept that rounding but
must still reject a transposed (x,y)->(y,x) read. Both cases are built from the real Friday values."""
import math, unittest

def samples_agree(expected, actual):
    """The patched rule from prepare_pilot_inputs.py."""
    for e, a in zip(expected, actual):
        same = (e['value'] is None and a['value'] is None) or (
            e['value'] is not None and a['value'] is not None and
            math.isclose(e['value'], a['value'], rel_tol=1e-12, abs_tol=1e-15))
        if e['fits_x_1based'] != a['fits_x_1based'] or e['fits_y_1based'] != a['fits_y_1based'] or not same:
            return False
    return True

PY_SAMPLES = [dict(fits_x_1based=2, fits_y_1based=5, value=0.012884007766842842),
              dict(fits_x_1based=5, fits_y_1based=2, value=-0.00876860786229372)]
R_SAMPLES_15DIGIT = [dict(fits_x_1based=2, fits_y_1based=5, value=0.0128840077668428),
                     dict(fits_x_1based=5, fits_y_1based=2, value=-0.00876860786229372)]
R_SAMPLES_TRANSPOSED = [dict(fits_x_1based=2, fits_y_1based=5, value=-0.00876860786229372),
                        dict(fits_x_1based=5, fits_y_1based=2, value=0.0128840077668428)]

class T(unittest.TestCase):
    def test_15_digit_rounding_accepted(self): self.assertTrue(samples_agree(PY_SAMPLES, R_SAMPLES_15DIGIT))
    def test_transposed_read_rejected(self): self.assertFalse(samples_agree(PY_SAMPLES, R_SAMPLES_TRANSPOSED))
    def test_exact_equality_would_have_failed(self): self.assertNotEqual(PY_SAMPLES[0]['value'], R_SAMPLES_15DIGIT[0]['value'])

if __name__ == '__main__': unittest.main(verbosity=2)
