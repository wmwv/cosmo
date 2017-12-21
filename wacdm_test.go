package cosmo

import (
	"gonum.org/v1/gonum/floats"
	"testing"
)

func TestWACDMCosmologyInterface(t *testing.T) {
	age_distance := func(cos FLRW) {
		z := 0.5
		age := cos.Age(z)
		dc := cos.ComovingDistance(z)
		_, _ = age, dc
	}

	cos := WACDM{Om0: 0.27, Ol0: 0.73, W0: -1, H0: 70, Tcmb0: 0.}
	age_distance(cos)
}

// TestE* tests that basic calculation of E
//   https://github.com/astropy/astropy/blob/master/astropy/cosmology/tests/test_cosmology.py
func TestWACDMELcdm(t *testing.T) {
	var exp, obs, tol float64
	cos := WACDM{Om0: 0.27, Ol0: 0.73, W0: -1, WA: 0, H0: 70, Tcmb0: 0.}

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
		t.Errorf("Failed flat waCDM E(z) test.  Expected %f, return %f",
			exp, obs)
	}

	exp = 1 / 1.7
	obs = cos.Einv(1.0)
	if !floats.EqualWithinAbs(obs, exp, tol) {
		t.Errorf("Failed flat waCDM Einv(z) test.  Expected %f, return %f",
			exp, obs)
	}
}

func TestWACDMDistanceModulus(t *testing.T) {
	var z_vec, exp_vec []float64
	var obs, tol float64
	cos := WACDM{Om0: 0.3, Ol0: 0.7, W0: -1.2, WA: 3, H0: 70, Tcmb0: 0.}

	tol = 1e-8
	// Calculated via
	//   from astropy.cosmology import w0waCDM
	//   z = np.asarray([0.5, 1.0, 2.0, 3.0])
	//   w0waCDM(70, 0.3, 0.7, -1.2, 3).distmod(z)
	z_vec = []float64{0.5, 1.0, 2.0, 3.0}
	exp_vec = []float64{42.20567831, 43.92122272, 45.57180818, 46.47483095}
	for i, z := range z_vec {
		obs = cos.DistanceModulus(z)
		if !floats.EqualWithinAbs(obs, exp_vec[i], tol) {
			t.Errorf("Failed flat waCDM luminosity distance test."+
				"  Expected %f, return %f",
				exp_vec[i], obs)
		}

	}
}

func TestWACDMLcdmDistanceModulus(t *testing.T) {
	var z_vec, exp_vec []float64
	var obs, tol float64
	cos := WACDM{Om0: 0.3, Ol0: 0.7, W0: -1.2, WA: 0, H0: 70, Tcmb0: 0.}

	tol = 1e-8
	// Calculated via
	//   from astropy.cosmology import w0waCDM
	//   z = np.asarray([0.5, 1.0, 2.0, 3.0])
	//   w0waCDM(70, 0.3, 0.7, -1.2, 0).distmod(z)
	z_vec = []float64{0.5, 1.0, 2.0, 3.0}
	exp_vec = []float64{42.32710911, 44.17957201, 46.03118144, 47.09228735}
	for i, z := range z_vec {
		obs = cos.DistanceModulus(z)
		if !floats.EqualWithinAbs(obs, exp_vec[i], tol) {
			t.Errorf("Failed flat waCDM luminosity distance test."+
				"  Expected %f, return %f",
				exp_vec[i], obs)
		}

	}
}

func TestWACDMLuminosityDistance(t *testing.T) {
	var z_vec, exp_vec []float64
	var obs, tol float64
	cos := WACDM{Om0: 0.3, Ol0: 0.7, W0: -1, WA: 2, H0: 70, Tcmb0: 0.}

	tol = 1e-6
	// Calculated via
	//   from astropy.cosmology import w0waCDM
	//   z = np.asarray([0.5, 1.0, 2.0, 3.0])
	//   w0waCDM(70, 0.3, 0.7, -1, 2).luminosity_distance(z)
	z_vec = []float64{0.5, 1.0, 2.0, 3.0}
	exp_vec = []float64{2723.8564247, 6067.27548444, 13407.53496469, 21024.89761811}
	for i, z := range z_vec {
		obs = cos.LuminosityDistance(z)
		if !floats.EqualWithinAbs(obs, exp_vec[i], tol) {
			t.Errorf("Failed flat waCDM luminosity distance test."+
				"  Expected %f, return %f",
				exp_vec[i], obs)
		}

	}
}

func TestWACDMAngularDiameterDistance(t *testing.T) {
	var z_vec, exp_vec []float64
	var obs, tol float64
	cos := WACDM{Om0: 0.3, Ol0: 0.7, W0: -0.8, WA: 2.5, H0: 70, Tcmb0: 0.}

	tol = 1e-6
	// Calculated via
	//   from astropy.cosmology import w0waCDM
	//   z = np.asarray([0.5, 1.0, 2.0, 3.0])
	//   w0waCDM(70, 0.3, 0.7, -0.8, 2.5).angular_diameter_distance(z)
	z_vec = []float64{0.5, 1.0, 2.0, 3.0}
	exp_vec = []float64{1155.52181127, 1393.61319898, 1282.08090454, 1073.63096224}
	for i, z := range z_vec {
		obs = cos.AngularDiameterDistance(z)
		if !floats.EqualWithinAbs(obs, exp_vec[i], tol) {
			t.Errorf("Failed flat waCDM angular diameter distance test."+
				"  Expected %f, return %f",
				exp_vec[i], obs)
		}

	}
}

func TestWACDMComovingTransverseDistance(t *testing.T) {
	var z_vec, exp_vec []float64
	var obs, tol float64
	cos := WACDM{Om0: 0.3, Ol0: 0.7, W0: -1.2, WA: -1.2, H0: 70, Tcmb0: 0.}

	tol = 1e-6
	// Calculated via
	//   from astropy.cosmology import w0waCDM
	//   z = np.asarray([0.5, 1.0, 2.0, 3.0])
	//   w0waCDM(70, 0.3, 0.7, -1.2, -1.2).comoving_transverse_distance(z)
	z_vec = []float64{0.5, 1.0, 2.0, 3.0}
	exp_vec = []float64{1985.54631561, 3533.91345688, 5524.66720808, 6731.56420461}
	for i, z := range z_vec {
		obs = cos.ComovingTransverseDistance(z)
		if !floats.EqualWithinAbs(obs, exp_vec[i], tol) {
			t.Errorf("Failed flat waCDM comoving transverse distance test."+
				"  Expected %f, return %f",
				exp_vec[i], obs)
		}

	}
}

func TestWACDMComovingDistanceZ1Z2Integrate(t *testing.T) {
	var z_vec, exp_vec []float64
	var obs, tol float64
	cos := WACDM{Om0: 0.3, Ol0: 0.7, W0: -0.9, WA: 3.5, H0: 70, Tcmb0: 0.}

	tol = 1e-6
	// Calculated via
	//   from astropy.cosmology import w0waCDM
	//   z = np.asarray([0.5, 1.0, 2.0, 3.0])
	//   w0waCDM(70, 0.3, 0.7, -0.9, 3.5)._comoving_distance_z1z2(0, z)
	z_vec = []float64{0.5, 1.0, 2.0, 3.0}
	exp_vec = []float64{1726.71519955, 2709.17698433, 3538.63486291, 3798.28908226}
	for i, z := range z_vec {
		obs = cos.ComovingDistanceZ1Z2Integrate(0, z)
		if !floats.EqualWithinAbs(obs, exp_vec[i], tol) {
			t.Errorf("Failed flat waCDM comoving distance elliptic test."+
				"  Expected %f, return %f",
				exp_vec[i], obs)
		}

	}
}

func TestWACDMLookbackTime(t *testing.T) {
	var z_vec, exp_vec []float64
	var obs, tol float64
	cos := WACDM{Om0: 0.3, Ol0: 0.7, W0: -0.9, WA: 3.5, H0: 70, Tcmb0: 0.}

	tol = 1e-6
	z_vec = []float64{0.5, 1.0, 2.0, 3.0}
	// Calculated via
	//   from astropy.cosmology import w0waCDM
	//   z = np.asarray([0.5, 1.0, 2.0, 3.0])
	//   w0waCDM(70, 0.3, 0.7, -0.9, 3.5).lookback_time(z)
	exp_vec = []float64{4.64427098, 6.51439755, 7.65559243, 7.90553458}
	for i, z := range z_vec {
		obs = cos.LookbackTime(z)
		if !floats.EqualWithinAbs(obs, exp_vec[i], tol) {
			t.Errorf("Failed flat waCDM w0=%.1f lookback time test."+
				"  Expected %f, return %f",
				cos.W0, exp_vec[i], obs)
		}
	}
}

func TestWACDMLookbackTimeIntegrate(t *testing.T) {
	var z_vec, exp_vec []float64
	var obs, tol float64
	cos := WACDM{Om0: 0.3, Ol0: 0.7, W0: -1.1, WA: 2.8, H0: 70, Tcmb0: 0.}

	tol = 1e-6
	// Calculated via
	//   from astropy.cosmology import w0waCDM
	//   z = np.asarray([0.5, 1.0, 2.0, 3.0])
	//   w0waCDM(70, 0.3, 0.7, -1.1, 2.8).lookback_time(z)
	z_vec = []float64{0.5, 1.0, 2.0, 3.0}
	exp_vec = []float64{4.86957219, 7.09750717, 8.81344126, 9.36817438}
	for i, z := range z_vec {
		obs = cos.LookbackTimeIntegrate(z)
		if !floats.EqualWithinAbs(obs, exp_vec[i], tol) {
			t.Errorf("Failed flat waCDM w0=%.1f lookback time integrate test."+
				"  Expected %f, return %f",
				cos.W0, exp_vec[i], obs)
		}
	}
}

func TestWACDMLookbackTimeOM(t *testing.T) {
	var z_vec, exp_vec []float64
	var obs, tol float64
	cos := WACDM{Om0: 0.3, Ol0: 0., W0: -0.9, WA: 2, H0: 70, Tcmb0: 0.}

	tol = 1e-6
	// Calculated via
	//   from astropy.cosmology import w0waCDM
	//   z = np.asarray([0.5, 1.0, 2.0, 3.0])
	//   w0waCDM(70, 0.3, 0., -0.9, 2).lookback_time(z)
	z_vec = []float64{0.5, 1.0, 2.0, 3.0}
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

func TestWACDMLookbackTimeOL(t *testing.T) {
	var z_vec, exp_vec []float64
	var obs, tol float64
	cos := WACDM{Om0: 0., Ol0: 0.5, W0: -1, WA: 0, H0: 70, Tcmb0: 0.}

	tol = 1e-6
	// Calculated via
	//   from astropy.cosmology import w0waCDM
	//   z = np.asarray([0.5, 1.0, 2.0, 3.0])
	//   w0waCDM(70, 0., 0.5, -1, 0).lookback_time(z)
	z_vec = []float64{0.5, 1.0, 2.0, 3.0}
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

func TestWACDMAge(t *testing.T) {
	var z_vec, exp_vec []float64
	var obs, tol float64
	cos := WACDM{Om0: 0.3, Ol0: 0.6, W0: -0.6, WA: 3.5, H0: 70, Tcmb0: 0.}

	tol = 1e-6

	// Calculated via
	//   import astropy.cosmology
	//   z = [0.5, 1.0, 2.0, 3.0]
	//   astropy.cosmology.w0waCDM(70, 0.3, 0.6, -0.6, 3.5).age(z)
	z_vec = []float64{0.5, 1.0, 2.0, 3.0}
	exp_vec = []float64{2.70980463, 1.08619498, 0.21688951, 0.058307}

	for i, z := range z_vec {
		obs = cos.Age(z)
		if !floats.EqualWithinAbs(obs, exp_vec[i], tol) {
			t.Errorf("Failed (OM, OL) = (%.2f, %.2f) LCDM age test."+
				"  Expected %f, return %f",
				cos.Om0, cos.Ol0, exp_vec[i], obs)
		}
	}
}
