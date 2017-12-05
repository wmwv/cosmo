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
	cos := Cosmology{Om0: 0.27, Ol0: 0.73, Ok0: 0., H0: 70, W0: -1.0, Tcmb0: 0.}

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
	cos := Cosmology{Om0: 0.3, Ol0: 0.7, Ok0: 0., H0: 70, W0: -1.0, Tcmb0: 0.}

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
	cos := Cosmology{Om0: 0.3, Ol0: 0.7, Ok0: 0., H0: 70, W0: -1.0, Tcmb0: 0.}

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

func TestAngularDiameterDistance(t *testing.T) {
	var z_vec, exp_vec []float64
	var obs, tol float64
	cos := Cosmology{Om0: 0.3, Ol0: 0.7, Ok0: 0., H0: 70, W0: -1.0, Tcmb0: 0.}

	tol = 1e-6
	//  z_vec = []float64{0.2, 0.4, 0.9, 1.2}
	//  exp_vec = []float64{971.667, 2141.67, 5685.96, 8107.41}
	z_vec = []float64{0.5, 1.0, 2.0, 3.0}
	exp_vec = []float64{1259.08359729, 1651.91440294, 1726.62069147, 1588.92135907}
	for i, z := range z_vec {
		obs = cos.AngularDiameterDistance(z)
		if !floats.EqualWithinAbs(obs, exp_vec[i], tol) {
			t.Errorf("Failed flat LCDM angular diameter distance test."+
				"  Expected %f, return %f",
				exp_vec[i], obs)
		}

	}
}

func TestComovingTransverseDistance(t *testing.T) {
	var z_vec, exp_vec []float64
	var obs, tol float64
	cos := Cosmology{Om0: 0.3, Ol0: 0.7, Ok0: 0., H0: 70, W0: -1.0, Tcmb0: 0.}

	tol = 1e-6
	//  z_vec = []float64{0.2, 0.4, 0.9, 1.2}
	//  exp_vec = []float64{971.667, 2141.67, 5685.96, 8107.41}
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

func TestComovingDistanceZ1Z2Integrate(t *testing.T) {
	var z_vec, exp_vec []float64
	var obs, tol float64
	cos := Cosmology{Om0: 0.3, Ol0: 0.7, Ok0: 0., H0: 70, W0: -1.0, Tcmb0: 0.}

	tol = 1e-6
	//  z_vec = []float64{0.2, 0.4, 0.9, 1.2}
	//  exp_vec = []float64{971.667, 2141.67, 5685.96, 8107.41}
	z_vec = []float64{0.5, 1.0, 2.0, 3.0}
	exp_vec = []float64{1888.62539593, 3303.82880589, 5179.86207441, 6355.6854363}
	for i, z := range z_vec {
		obs = cos.ComovingDistanceZ1Z2Integrate(0, z)
		if !floats.EqualWithinAbs(obs, exp_vec[i], tol) {
			t.Errorf("Failed flat LCDM comoving distance elliptic test."+
				"  Expected %f, return %f",
				exp_vec[i], obs)
		}

	}
}

func TestComovingDistanceZ1Z2Elliptic(t *testing.T) {
	var z_vec, exp_vec []float64
	var obs, tol float64
	cos := Cosmology{Om0: 0.3, Ol0: 0.7, Ok0: 0., H0: 70, W0: -1.0, Tcmb0: 0.}

	tol = 1e-6
	//	z_vec = []float64{0.2, 0.4, 0.9, 1.2}
	//	exp_vec = []float64{971.667, 2141.67, 5685.96, 8107.41}
	z_vec = []float64{0.5, 1.0, 2.0, 3.0}
	exp_vec = []float64{1888.62539593, 3303.82880589, 5179.86207441, 6355.6854363}
	for i, z := range z_vec {
		obs = cos.ComovingDistanceZ1Z2Elliptic(0, z)
		if !floats.EqualWithinAbs(obs, exp_vec[i], tol) {
			t.Errorf("Failed flat LCDM comoving distance elliptic test."+
				"  Expected %f, return %f",
				exp_vec[i], obs)
		}

	}
}

func TestLookbackTime(t *testing.T) {
	var z_vec, exp_vec []float64
	var obs, tol float64
	cos := Cosmology{Om0: 0.3, Ol0: 0.7, Ok0: 0., H0: 70, W0: -1.0, Tcmb0: 0.}

	tol = 1e-6
	z_vec = []float64{0.5, 1.0, 2.0, 3.0}
	// Calculated via astropy.cosmology.FlatLambdaCDM(70, 0.3).lookback_time
	exp_vec = []float64{5.04063793, 7.715337, 10.24035689, 11.35445676}
	for i, z := range z_vec {
		obs = cos.LookbackTime(z)
		if !floats.EqualWithinAbs(obs, exp_vec[i], tol) {
			t.Errorf("Failed flat LCDM lookback time test."+
				"  Expected %f, return %f",
				exp_vec[i], obs)
		}
	}
}

func TestLookbackTimeIntegrate(t *testing.T) {
	var z_vec, exp_vec []float64
	var obs, tol float64
	cos := Cosmology{Om0: 0.3, Ol0: 0.7, Ok0: 0., H0: 70, W0: -1.0, Tcmb0: 0.}

	tol = 1e-6
	z_vec = []float64{0.5, 1.0, 2.0, 3.0}
	// Calculated via astropy.cosmology.FlatLambdaCDM(70, 0.3).lookback_time
	exp_vec = []float64{5.04063793, 7.715337, 10.24035689, 11.35445676}
	for i, z := range z_vec {
		obs = cos.LookbackTimeIntegrate(z)
		if !floats.EqualWithinAbs(obs, exp_vec[i], tol) {
			t.Errorf("Failed flat LCDM lookback time integrate test."+
				"  Expected %f, return %f",
				exp_vec[i], obs)
		}
	}
}

func TestLookbackTimeOM(t *testing.T) {
	var z_vec, exp_vec []float64
	var obs, tol float64
	cos := Cosmology{Om0: 0.3, Ol0: 0., Ok0: 0., H0: 70, W0: -1.0, Tcmb0: 0.}

	tol = 1e-6
	z_vec = []float64{0.5, 1.0, 2.0, 3.0}
	// Calculated via astropy.cosmology.FlatLambdaCDM(70, 0.3).lookback_time
	exp_vec = []float64{4.51471693, 6.62532254, 8.57486509, 9.45923582}
	for i, z := range z_vec {
		obs = cos.LookbackTimeOM(z)
		if !floats.EqualWithinAbs(obs, exp_vec[i], tol) {
			t.Errorf("Failed (OM, OL)=(0.3, 0) lookback time test."+
				"  Expected %f, return %f",
				exp_vec[i], obs)
		}
	}
}

func TestLookbackTimeOL(t *testing.T) {
	var z_vec, exp_vec []float64
	var obs, tol float64
	cos := Cosmology{Om0: 0., Ol0: 0.5, Ok0: 0., H0: 70, W0: -1.0, Tcmb0: 0.}

	tol = 1e-6
	z_vec = []float64{0.5, 1.0, 2.0, 3.0}
	// Calculated via astropy.cosmology.FlatLambdaCDM(70, 0.3).lookback_time
	exp_vec = []float64{5.0616361, 7.90494991, 10.94241739, 12.52244605}
	for i, z := range z_vec {
		obs = cos.LookbackTimeOL(z)
		if !floats.EqualWithinAbs(obs, exp_vec[i], tol) {
			t.Errorf("Failed (OM, OL)=(0, 0.5) lookback time test."+
				"  Expected %f, return %f",
				exp_vec[i], obs)
		}
	}
}

func TestAge(t *testing.T) {
	var z_vec, exp_vec []float64
	var obs, tol float64
	cos := Cosmology{Om0: 0.3, Ol0: 0.6, Ok0: 0., H0: 70, W0: -1.0, Tcmb0: 0.}

	tol = 1e-6

	z_vec = []float64{0.5, 1.0, 2.0, 3.0}
	// Calculated via
	//   import astropy.cosmology
	//   z = [0.5, 1.0, 2.0, 3.0]
	//   astropy.cosmology.LambdaCDM(70, 0.3, 0.6).age(z)
	exp_vec = []float64{8.11137578, 5.54558439, 3.13456008, 2.06445301}

	for i, z := range z_vec {
		obs = cos.Age(z)
		if !floats.EqualWithinAbs(obs, exp_vec[i], tol) {
			t.Errorf("Failed (OM, OL) = (%.2f, %.2f) LCDM age test."+
				"  Expected %f, return %f",
				cos.Om0, cos.Ol0, exp_vec[i], obs)
		}
	}
}

func TestAgeFlatLCDM(t *testing.T) {
	var z_vec, exp_vec []float64
	var obs, tol float64
	cos := Cosmology{Om0: 0.3, Ol0: 0.7, Ok0: 0., H0: 70, W0: -1.0, Tcmb0: 0.}

	tol = 1e-6
	z_vec = []float64{0.5, 1.0, 2.0, 3.0}

	// Calculated via
	//   import astropy.cosmology
	//   z = [0.5, 1.0, 2.0, 3.0]
	//   astropy.cosmology.FlatLambdaCDM(70, 0.3).age(z)
	exp_vec = []float64{8.42634602, 5.75164694, 3.22662706, 2.11252719}

	for i, z := range z_vec {
		obs = cos.Age(z)
		if !floats.EqualWithinAbs(obs, exp_vec[i], tol) {
			t.Errorf("Failed flat LCDM age test."+
				"  Expected %f, return %f",
				exp_vec[i], obs)
		}
	}
}

func TestAgeIntegrate(t *testing.T) {
	var z_vec, exp_vec []float64
	var obs, tol float64
	cos := Cosmology{Om0: 0.3, Ol0: 0.7, Ok0: 0., H0: 70, W0: -1.0, Tcmb0: 0.}

	tol = 1e-6
	z_vec = []float64{0.5, 1.0, 2.0, 3.0}

	// Calculated via
	//   import astropy.cosmology
	//   z = [0.5, 1.0, 2.0, 3.0]
	//   astropy.cosmology.FlatLambdaCDM(70, 0.3).age(z)
	exp_vec = []float64{8.42634602, 5.75164694, 3.22662706, 2.11252719}

	for i, z := range z_vec {
		obs = cos.AgeIntegrate(z)
		if !floats.EqualWithinAbs(obs, exp_vec[i], tol) {
			t.Errorf("Failed flat LCDM ageIntegrate test."+
				"  Expected %f, return %f",
				exp_vec[i], obs)
		}
	}
}

func TestAgeOM(t *testing.T) {
	var z_vec, exp_vec []float64
	var obs, tol float64
	var ageIntegrate float64
	cos := Cosmology{Om0: 0.3, Ol0: 0., Ok0: 0., H0: 70, W0: -1.0, Tcmb0: 0.}

	tol = 1e-6
	z_vec = []float64{0.5, 1.0, 2.0, 3.0}
	// Calculated via
	//   import astropy.cosmology
	//   z = [0.5, 1.0, 2.0, 3.0]
	//   astropy.cosmology.LambdaCDM(70, 0.3, 0.).age(z)
	exp_vec = []float64{6.78287955, 4.67227393, 2.72273139, 1.83836065}
	for i, z := range z_vec {
		obs = cos.AgeOM(z)
		if !floats.EqualWithinAbs(obs, exp_vec[i], tol) {
			t.Errorf("Failed (OM, OL)=(0.3, 0) age test."+
				"  Expected %f, return %f",
				exp_vec[i], obs)
		}
		ageIntegrate = cos.AgeIntegrate(z)
		if !floats.EqualWithinAbs(obs, ageIntegrate, tol) {
			t.Errorf("Failed (OM, OL)=(0.3, 0) comparison "+
				"between AgeOM (%f) and AgeIntegrate (%f).",
				obs, ageIntegrate)
		}
	}
}

func TestAgeOL(t *testing.T) {
	var z_vec, exp_vec []float64
	var obs, tol float64
	cos := Cosmology{Om0: 0., Ol0: 0.5, Ok0: 0., H0: 70, W0: -1.0, Tcmb0: 0.}

	tol = 1e-6
	z_vec = []float64{0.5, 1.0, 2.0, 3.0}
	// Calculated via astropy.cosmology.FlatLambdaCDM(70, 0.3).lookback_time
	exp_vec = []float64{12.34935796, 9.50604415, 6.46857667, 4.88854801}
	for i, z := range z_vec {
		obs = cos.AgeOL(z)
		if !floats.EqualWithinAbs(obs, exp_vec[i], tol) {
			t.Errorf("Failed (OM, OL)=(0, 0.5) age test."+
				"  Expected %f, return %f",
				exp_vec[i], obs)
		}
	}
}

// Analytic case of Omega_Lambda = 0
func TestEOm(t *testing.T) {
	var z_vec []float64
	var obs, exp, tol float64
	cos := Cosmology{Om0: 1.0, Ol0: 0., Ok0: 0., H0: 70, W0: -1.0, Tcmb0: 0.}
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
	cos := Cosmology{Om0: 0.27, Ol0: 0.73, Ok0: 0., H0: 70, W0: -1.0, Tcmb0: 0.}
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
	cos := Cosmology{Om0: 0.27, Ol0: 0.73, Ok0: 0., H0: 70, W0: -0.9, Tcmb0: 0.}

	// Check value of E(z=1.0)
	exp := 1.7489240754
	obs := cos.E(1.0)
	tol := 1e-9
	if !floats.EqualWithinAbs(obs, exp, tol) {
		t.Errorf("Expected %f, return %f", exp, obs)
	}
}
