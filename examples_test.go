package cosmo

import (
	"gonum.org/v1/gonum/floats"
	"testing"
)

// TestE* tests that basic calculation of E
//   https://github.com/astropy/astropy/blob/master/astropy/cosmology/tests/test_cosmology.py
func TestELcdm(t *testing.T) {
	var exp, obs, tol float64
	cos := Cosmology{Om0: 0.27, Ol0: 0.73, Ok0: 0., H0: 70, w0: -1.0, Tcmb0: 0.}

	// Check value of E(z=1.0)
	//   OM, OL, OK, z = 0.27, 0.73, 0.0, 1.0
	//   sqrt(OM*(1+z)**3 + OK * (1+z)**2 + OL)
	//   sqrt(0.27*(1+1.0)**3 + 0.0 * (1+1.0)**2 + 0.73)
	//   sqrt(0.27*8 + 0 + 0.73)
	//   sqrt(2.89)
	exp = 1.7
	obs = cos.E(1.0)
	tol = 1e-9
	if !floats.EqualWithinAbs(obs, exp, tol) {
		t.Errorf("Failed flat LCDM E(z) test.  Expected %f, return %f",
			exp, obs)
	}

	exp = 1 / 1.7
	obs = cos.Einv(1.0)
	if !floats.EqualWithinAbs(obs, exp, tol) {
		t.Errorf("Failed flat LCDM Einv(z) test.  Expected %f, return %f",
			exp, obs)
	}

}

func TestEvecLcdm(t *testing.T) {
	cos := Cosmology{Om0: 0.27, Ol0: 0.73, Ok0: 0., H0: 70, w0: -1.0, Tcmb0: 0.}
	// Check array
	z := []float64{0.5, 1.0}
	// Calculated using astropy.cosmology.FlatLambdaCDM (v1.3.2)
	exp := []float64{1.2811127975318957, 1.7}
	tol := 1e-9
	obs := cos.Evec(z)
	if !floats.EqualApprox(obs, exp, tol) {
		t.Errorf("Failed array float LCDM test.  Expected %v, return %v", exp, obs)
	}

}

// Check value of E(z = [0.5, 1.0])
// Testing both value and array
func TestEwCDM(t *testing.T) {
	cos := Cosmology{Om0: 0.27, Ol0: 0.73, Ok0: 0., H0: 70, w0: -0.9, Tcmb0: 0.}

	// Check value of E(z=1.0)
	exp := 1.7489240754
	obs := cos.E(1.0)
	tol := 1e-9
	if !floats.EqualWithinAbs(obs, exp, tol) {
		t.Errorf("Expected %f, return %f", exp, obs)
	}
}
