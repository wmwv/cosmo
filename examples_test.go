package cosmo

import (
	"gonum.org/v1/gonum/floats"
	"testing"
)

// TestE tests that basic calculation of E
//   https://github.com/astropy/astropy/blob/master/astropy/cosmology/tests/test_cosmology.py
func TestE(t *testing.T) {
	cos := Cosmology{Om0: 0.27, Ol0: 0.73, Ok0: 0., H0: 70, w0: -0.9, Tcmb0: 0.}

	// Check value of E(z=1.0)
	exp := 1.7489240754
	obs := cos.E(1.0)
	tol := 1e-9
	if !floats.EqualWithinAbs(obs, exp, tol) {
		t.Errorf("Expected %f, return %f", exp, obs)
	}
}
