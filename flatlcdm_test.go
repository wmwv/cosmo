package cosmo

import (
	"gonum.org/v1/gonum/floats"
	"math"
	"testing"
)

func TestFlatLCDMCosmologyInterface(t *testing.T) {
	age_distance := func(cos FLRW) {
		z := 0.5
		age := cos.Age(z)
		dc := cos.ComovingDistance(z)
		_, _ = age, dc
	}

	cos := FlatLCDM{Om0: 0.27, H0: 70, Tcmb0: 0.}
	age_distance(cos)
}

// TestE* tests that basic calculation of E
//   https://github.com/astropy/astropy/blob/master/astropy/cosmology/tests/test_cosmology.py
func TestFlatLCDME(t *testing.T) {
	var exp, obs, tol float64
	cos := FlatLCDM{Om0: 0.27, H0: 70, Tcmb0: 0.}

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

func TestFlatLCDMDistanceModulus(t *testing.T) {
	var z_vec, exp_vec []float64
	var obs, tol float64
	cos := FlatLCDM{Om0: 0.3, H0: 70, Tcmb0: 0.}

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

func TestFlatLCDMLuminosityDistance(t *testing.T) {
	var z_vec, exp_vec []float64
	var obs, tol float64
	cos := FlatLCDM{Om0: 0.3, H0: 70, Tcmb0: 0.}

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

func TestFlatLCDMAngularDiameterDistance(t *testing.T) {
	var z_vec, exp_vec []float64
	var obs, tol float64
	cos := FlatLCDM{Om0: 0.3, H0: 70, Tcmb0: 0.}

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

func TestFlatLCDMComovingTransverseDistance(t *testing.T) {
	var z_vec, exp_vec []float64
	var obs, tol float64
	cos := FlatLCDM{Om0: 0.3, H0: 70, Tcmb0: 0.}

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

func TestFlatLCDMComovingDistanceZ1Z2Integrate(t *testing.T) {
	var z_vec, exp_vec []float64
	var obs, tol float64
	cos := FlatLCDM{Om0: 0.3, H0: 70, Tcmb0: 0.}

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

func TestFlatLCDMComovingDistanceZ1Z2Elliptic(t *testing.T) {
	var z_vec, exp_vec []float64
	var obs, tol float64
	cos := FlatLCDM{Om0: 0.3, H0: 70, Tcmb0: 0.}

	tol = 1e-6
	//  z_vec = []float64{0.2, 0.4, 0.9, 1.2}
	//  exp_vec = []float64{971.667, 2141.67, 5685.96, 8107.41}
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

func TestFlatLCDMLookbackTime(t *testing.T) {
	var z_vec, exp_vec []float64
	var obs, tol float64
	cos := FlatLCDM{Om0: 0.3, H0: 70, Tcmb0: 0.}

	tol = 1e-6
	z_vec = []float64{0.5, 1.0, 2.0, 3.0}
	// Calculated via astropy.cosmology.FlatFlatLCDM(70, 0.3).lookback_time
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

func TestFlatLCDMLookbackTimeIntegrate(t *testing.T) {
	var z_vec, exp_vec []float64
	var obs, tol float64
	cos := FlatLCDM{Om0: 0.3, H0: 70, Tcmb0: 0.}

	tol = 1e-6
	z_vec = []float64{0.5, 1.0, 2.0, 3.0}
	// Calculated via astropy.cosmology.FlatFlatLCDM(70, 0.3).lookback_time
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

func TestFlatLCDMAgeFlatLCDM(t *testing.T) {
	var z_vec, exp_vec []float64
	var obs, tol float64
	cos := FlatLCDM{Om0: 0.3, H0: 70, Tcmb0: 0.}

	tol = 1e-6
	z_vec = []float64{0.5, 1.0, 2.0, 3.0}

	// Calculated via
	//   import astropy.cosmology
	//   z = [0.5, 1.0, 2.0, 3.0]
	//   astropy.cosmology.FlatFlatLCDM(70, 0.3).age(z)
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

func TestFlatLCDMAgeIntegrate(t *testing.T) {
	var z_vec, exp_vec []float64
	var obs, tol float64
	cos := FlatLCDM{Om0: 0.3, H0: 70, Tcmb0: 0.}

	tol = 1e-6
	z_vec = []float64{0.5, 1.0, 2.0, 3.0}

	// Calculated via
	//   import astropy.cosmology
	//   z = [0.5, 1.0, 2.0, 3.0]
	//   astropy.cosmology.FlatFlatLCDM(70, 0.3).age(z)
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

// Analytic case of Omega_Lambda = 0
func TestFlatLCDMEOm(t *testing.T) {
	var z_vec []float64
	var obs, exp, tol float64
	cos := FlatLCDM{Om0: 1.0, H0: 70, Tcmb0: 0.}
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

func TestFlatLCDMEvecLcdm(t *testing.T) {
	cos := FlatLCDM{Om0: 0.27, H0: 70, Tcmb0: 0.}
	// Check array
	z := []float64{0.5, 1.0}
	// Calculated using astropy.cosmology.FlatFlatLCDM (v1.3.2)
	exp := []float64{1.2811127975318957, 1.7}
	tol := 1e-9
	obs := cos.Evec(z)
	if !floats.EqualApprox(obs, exp, tol) {
		t.Errorf("Failed array float LCDM test.  Expected %v, return %v", exp, obs)
	}

}
