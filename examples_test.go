package cosmo

import (
	"gonum.org/v1/gonum/floats"
	"math"
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

func TestDistanceModulus(t *testing.T) {
	var z_vec, exp_vec []float64
	var obs, tol float64
	cos := Cosmology{Om0: 0.3, Ol0: 0.7, Ok0: 0., H0: 70, w0: -1.0, Tcmb0: 0.}

	tol = 1e-8
	//  z_vec = []float64{0.2, 0.4, 0.9, 1.2}
	//  exp_vec = []float64{971.667, 2141.67, 5685.96, 8107.41}
	z_vec = []float64{0.5, 1.0, 2.0, 3.0}
	exp_vec = []float64{42.26118542, 44.10023766, 45.95719725, 47.02611193}
	for i, z := range z_vec {
		obs = cos.DistanceModulus(z)
		if !floats.EqualWithinAbs(obs, exp_vec[i], tol) {
			t.Errorf("Failed flat LCDM luminosity distance test."+
				"  Expected %f, return %f",
				exp_vec[i], obs)
		}

	}
}

func TestLuminosityDistance(t *testing.T) {
	var z_vec, exp_vec []float64
	var obs, tol float64
	cos := Cosmology{Om0: 0.3, Ol0: 0.7, Ok0: 0., H0: 70, w0: -1.0, Tcmb0: 0.}

	tol = 1e-6
	//  z_vec = []float64{0.2, 0.4, 0.9, 1.2}
	//  exp_vec = []float64{971.667, 2141.67, 5685.96, 8107.41}
	z_vec = []float64{0.5, 1.0, 2.0, 3.0}
	exp_vec = []float64{2832.9380939, 6607.65761177, 15539.58622323, 25422.74174519}
	for i, z := range z_vec {
		obs = cos.LuminosityDistance(z)
		if !floats.EqualWithinAbs(obs, exp_vec[i], tol) {
			t.Errorf("Failed flat LCDM luminosity distance test."+
				"  Expected %f, return %f",
				exp_vec[i], obs)
		}

	}
}

func TestComovingTransverseDistance(t *testing.T) {
	var z_vec, exp_vec []float64
	var obs, tol float64
	cos := Cosmology{Om0: 0.3, Ol0: 0.7, Ok0: 0., H0: 70, w0: -1.0, Tcmb0: 0.}

	tol = 1e-6
	//	z_vec = []float64{0.2, 0.4, 0.9, 1.2}
	//	exp_vec = []float64{971.667, 2141.67, 5685.96, 8107.41}
	z_vec = []float64{0.5, 1.0, 2.0, 3.0}
	exp_vec = []float64{1888.62539593, 3303.82880589, 5179.86207441, 6355.6854363}
	for i, z := range z_vec {
		obs = cos.ComovingTransverseDistance(z)
		if !floats.EqualWithinAbs(obs, exp_vec[i], tol) {
			t.Errorf("Failed flat LCDM comoving transverse distance test."+
				"  Expected %f, return %f",
				exp_vec[i], obs)
		}

	}
}

// Analytic case of Omega_Lambda = 0
func TestEOm(t *testing.T) {
	var z_vec []float64
	var obs, exp, tol float64
	cos := Cosmology{Om0: 1.0, Ol0: 0., Ok0: 0., H0: 70, w0: -1.0, Tcmb0: 0.}
	tol = 1e-9
	z_vec = []float64{1.0, 10.0, 500.0, 1000.0}
	hubbleDistance := SpeedOfLightKmS / cos.H0
	for _, z := range z_vec {
		exp = 2.0 * hubbleDistance * (1 - math.Sqrt(1/(1+z)))
		obs = cos.ComovingDistance(z)
		if !floats.EqualWithinAbs(obs, exp, tol) {
			t.Errorf("Failed OM, OL = (1, 0) analytic comoving distance test."+
				"  Expected %f, return %f",
				exp, obs)
		}
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
