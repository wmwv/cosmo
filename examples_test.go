package cosmo

import (
	"gonum.org/v1/gonum/floats"
	"testing"
)

// TestE tests that basic calculation of E
//   https://github.com/astropy/astropy/blob/master/astropy/cosmology/tests/test_cosmology.py
func TestE(t *testing.T) {
	cos := Cosmology{Om0: 0.27, Ol0: 0.73, Ok0: 0., H0: 70, w0: -1.0, Tcmb0: 0.}

	// Check value of E(z=1.0)
    //   OM, OL, OK, z = 0.27, 0.73, 0.0, 1.0
    //   sqrt(OM*(1+z)**3 + OK * (1+z)**2 + OL)
    //   sqrt(0.27*(1+1.0)**3 + 0.0 * (1+1.0)**2 + 0.73)
    //   sqrt(0.27*8 + 0 + 0.73)
    //   sqrt(2.89)
	exp := 1.7
	obs := cos.E(1.0)
	tol := 1e-9
	if !floats.EqualWithinAbs(obs, exp, tol) {
		t.Errorf("Expected %f, return %f", exp, obs)
	}
}
