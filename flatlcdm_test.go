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
	var z, exp, tol float64
	cos := FlatLCDM{Om0: 0.27, H0: 70, Tcmb0: 0.}

	// Check value of E(z=1.0)
	//   OM, OL, OK, z = 0.27, 0.73, 0.0, 1.0
	//   sqrt(OM*(1+z)**3 + OK * (1+z)**2 + OL)
	//   sqrt(0.27*(1+1.0)**3 + 0.0 * (1+1.0)**2 + 0.73)
	//   sqrt(0.27*8 + 0 + 0.73)
	//   sqrt(2.89)
	z = 1.0
	exp = 1.7
	tol = 1e-9
	runTest(cos.E, z, exp, tol, t, 0)

	exp = 1 / 1.7
	runTest(cos.Einv, z, exp, tol, t, 0)
}

func TestFlatLCDMDistanceModulus(t *testing.T) {
	var z_vec, exp_vec []float64
	cos := FlatLCDM{Om0: 0.3, H0: 70, Tcmb0: 0.}

	tol := 1e-8
	//  z_vec = []float64{0.2, 0.4, 0.9, 1.2}
	//  exp_vec = []float64{971.667, 2141.67, 5685.96, 8107.41}
	z_vec = []float64{0.5, 1.0, 2.0, 3.0}
	exp_vec = []float64{42.26118542, 44.10023766, 45.95719725, 47.02611193}
	runTests(cos.DistanceModulus, z_vec, exp_vec, tol, t)
}

func TestFlatLCDMLuminosityDistance(t *testing.T) {
	var z_vec, exp_vec []float64
	cos := FlatLCDM{Om0: 0.3, H0: 70, Tcmb0: 0.}

	tol := 1e-6
	//  z_vec = []float64{0.2, 0.4, 0.9, 1.2}
	//  exp_vec = []float64{971.667, 2141.67, 5685.96, 8107.41}
	z_vec = []float64{0.5, 1.0, 2.0, 3.0}
	exp_vec = []float64{2832.9380939, 6607.65761177, 15539.58622323, 25422.74174519}
	runTests(cos.LuminosityDistance, z_vec, exp_vec, tol, t)
}

func TestFlatLCDMAngularDiameterDistance(t *testing.T) {
	var z_vec, exp_vec []float64
	cos := FlatLCDM{Om0: 0.3, H0: 70, Tcmb0: 0.}

	tol := 1e-6
	//  z_vec = []float64{0.2, 0.4, 0.9, 1.2}
	//  exp_vec = []float64{971.667, 2141.67, 5685.96, 8107.41}
	z_vec = []float64{0.5, 1.0, 2.0, 3.0}
	exp_vec = []float64{1259.08359729, 1651.91440294, 1726.62069147, 1588.92135907}
	runTests(cos.AngularDiameterDistance, z_vec, exp_vec, tol, t)
}

func TestFlatLCDMComovingTransverseDistance(t *testing.T) {
	var z_vec, exp_vec []float64
	cos := FlatLCDM{Om0: 0.3, H0: 70, Tcmb0: 0.}

	tol := 1e-6
	//  z_vec = []float64{0.2, 0.4, 0.9, 1.2}
	//  exp_vec = []float64{971.667, 2141.67, 5685.96, 8107.41}
	z_vec = []float64{0.5, 1.0, 2.0, 3.0}
	exp_vec = []float64{1888.62539593, 3303.82880589, 5179.86207441, 6355.6854363}
	runTests(cos.ComovingTransverseDistance, z_vec, exp_vec, tol, t)
}

func TestFlatLCDMComovingDistanceZ1Z2Integrate(t *testing.T) {
	var z_vec, exp_vec []float64
	cos := FlatLCDM{Om0: 0.3, H0: 70, Tcmb0: 0.}

	tol := 1e-6
	//  z_vec = []float64{0.2, 0.4, 0.9, 1.2}
	//  exp_vec = []float64{971.667, 2141.67, 5685.96, 8107.41}
	z_vec = []float64{0.5, 1.0, 2.0, 3.0}
	exp_vec = []float64{1888.62539593, 3303.82880589, 5179.86207441, 6355.6854363}
	runTestsZ0Z2(cos.ComovingDistanceZ1Z2Integrate, z_vec, exp_vec, tol, t)
}

func TestFlatLCDMComovingDistanceZ1Z2Elliptic(t *testing.T) {
	var z_vec, exp_vec []float64
	cos := FlatLCDM{Om0: 0.3, H0: 70, Tcmb0: 0.}

	tol := 1e-6
	//  z_vec = []float64{0.2, 0.4, 0.9, 1.2}
	//  exp_vec = []float64{971.667, 2141.67, 5685.96, 8107.41}
	z_vec = []float64{0.5, 1.0, 2.0, 3.0}
	exp_vec = []float64{1888.62539593, 3303.82880589, 5179.86207441, 6355.6854363}
	runTestsZ0Z2(cos.ComovingDistanceZ1Z2Elliptic, z_vec, exp_vec, tol, t)
}

func TestFlatLCDMLookbackTime(t *testing.T) {
	var z_vec, exp_vec []float64
	cos := FlatLCDM{Om0: 0.3, H0: 70, Tcmb0: 0.}

	tol := 1e-6
	z_vec = []float64{0.5, 1.0, 2.0, 3.0}
	// Calculated via astropy.cosmology.FlatFlatLCDM(70, 0.3).lookback_time
	exp_vec = []float64{5.04063793, 7.715337, 10.24035689, 11.35445676}
	runTests(cos.LookbackTime, z_vec, exp_vec, tol, t)
}

func TestFlatLCDMLookbackTimeIntegrate(t *testing.T) {
	var z_vec, exp_vec []float64
	cos := FlatLCDM{Om0: 0.3, H0: 70, Tcmb0: 0.}

	tol := 1e-6
	z_vec = []float64{0.5, 1.0, 2.0, 3.0}
	// Calculated via astropy.cosmology.FlatFlatLCDM(70, 0.3).lookback_time
	exp_vec = []float64{5.04063793, 7.715337, 10.24035689, 11.35445676}
	runTests(cos.LookbackTimeIntegrate, z_vec, exp_vec, tol, t)
}

func TestFlatLCDMAgeFlatLCDM(t *testing.T) {
	var z_vec, exp_vec []float64
	cos := FlatLCDM{Om0: 0.3, H0: 70, Tcmb0: 0.}

	tol := 1e-6
	z_vec = []float64{0.5, 1.0, 2.0, 3.0}

	// Calculated via
	//   import astropy.cosmology
	//   z = [0.5, 1.0, 2.0, 3.0]
	//   astropy.cosmology.FlatFlatLCDM(70, 0.3).age(z)
	exp_vec = []float64{8.42634602, 5.75164694, 3.22662706, 2.11252719}
	runTests(cos.Age, z_vec, exp_vec, tol, t)
}

func TestFlatLCDMAgeIntegrate(t *testing.T) {
	var z_vec, exp_vec []float64
	cos := FlatLCDM{Om0: 0.3, H0: 70, Tcmb0: 0.}

	tol := 1e-6
	z_vec = []float64{0.5, 1.0, 2.0, 3.0}

	// Calculated via
	//   import astropy.cosmology
	//   z = [0.5, 1.0, 2.0, 3.0]
	//   astropy.cosmology.FlatFlatLCDM(70, 0.3).age(z)
	exp_vec = []float64{8.42634602, 5.75164694, 3.22662706, 2.11252719}
	runTests(cos.AgeIntegrate, z_vec, exp_vec, tol, t)
}

// Analytic case of Omega_Lambda = 0
func TestFlatLCDMEOm(t *testing.T) {
	var z_vec, exp_vec []float64
	cos := FlatLCDM{Om0: 1.0, H0: 70, Tcmb0: 0.}
	tol := 1e-9
	z_vec = []float64{1.0, 10.0, 500.0, 1000.0}
	hubbleDistance := SpeedOfLightKmS / cos.H0
	exp_vec = make([]float64, len(z_vec))
	for i, z := range z_vec {
		exp_vec[i] = 2.0 * hubbleDistance * (1 - math.Sqrt(1/(1+z)))
	}
	runTests(cos.ComovingDistance, z_vec, exp_vec, tol, t)
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
