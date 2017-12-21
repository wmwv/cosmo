package cosmo

import (
	"gonum.org/v1/gonum/floats"
	"math"
	"testing"
)

var zLambdaCDM = []float64{0.5, 1.0, 2.0, 3.0}

// Calculated via Python AstroPy
//   from astropy.cosmology import LambdaCDM
//   z = np.asarray([0.5, 1.0, 2.0, 3.0])
var answers = map[string][]float64{}

func TestLambdaCDMCosmologyInterface(t *testing.T) {
	age_distance := func(cos FLRW) {
		z := 0.5
		age := cos.Age(z)
		dc := cos.ComovingDistance(z)
		_, _ = age, dc
	}

	cos := LambdaCDM{Om0: 0.27, Ol0: 0.73, H0: 70, Tcmb0: 0.}
	age_distance(cos)
}

// TestE* tests that basic calculation of E
//   https://github.com/astropy/astropy/blob/master/astropy/cosmology/tests/test_cosmology.py
func TestLambdaCDMELcdm(t *testing.T) {
	var z, exp, tol float64
	cos := LambdaCDM{Om0: 0.27, Ol0: 0.73, H0: 70, Tcmb0: 0.}

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

func TestLambdaCDMDistanceModulus(t *testing.T) {
	cos := LambdaCDM{Om0: 0.3, Ol0: 0.7, H0: 70, Tcmb0: 0.}
	exp_vec := []float64{42.26118542, 44.10023766, 45.95719725, 47.02611193}
	runTests(cos.DistanceModulus, zLambdaCDM, exp_vec, distTol, t)
}

func TestLambdaCDMLuminosityDistanceFlat(t *testing.T) {
	cos := LambdaCDM{Om0: 0.3, Ol0: 0.7, H0: 70, Tcmb0: 0.}
	//   LambdaCDM(70, 0.3, 0.7).luminosity_distance(z)
	exp_vec := []float64{2832.9380939, 6607.65761177, 15539.58622323, 25422.74174519}
	runTests(cos.LuminosityDistance, zLambdaCDM, exp_vec, distTol, t)
}

func TestLambdaCDMLuminosityDistanceNonflat(t *testing.T) {
	cos := LambdaCDM{Om0: 0.3, Ol0: 0.6, H0: 70, Tcmb0: 0.}
	//   LambdaCDM(70, 0.3, 0.6).luminosity_distance(z)
	exp_vec := []float64{2787.51504671, 6479.83450953, 15347.21516211, 25369.7240234}
	runTests(cos.LuminosityDistance, zLambdaCDM, exp_vec, distTol, t)
}

func TestLambdaCDMAngularDiameterDistance(t *testing.T) {
	cos := LambdaCDM{Om0: 0.3, Ol0: 0.7, H0: 70, Tcmb0: 0.}
	exp_vec := []float64{1259.08359729, 1651.91440294, 1726.62069147, 1588.92135907}
	runTests(cos.AngularDiameterDistance, zLambdaCDM, exp_vec, distTol, t)
}

func TestLambdaCDMComovingTransverseDistance(t *testing.T) {
	cos := LambdaCDM{Om0: 0.3, Ol0: 0.7, H0: 70, Tcmb0: 0.}
	exp_vec := []float64{1888.62539593, 3303.82880589, 5179.86207441, 6355.6854363}
	runTests(cos.ComovingTransverseDistance, zLambdaCDM, exp_vec, distTol, t)
}

func TestLambdaCDMComovingDistanceZ1Z2Integrate(t *testing.T) {
	cos := LambdaCDM{Om0: 0.3, Ol0: 0.7, H0: 70, Tcmb0: 0.}
	exp_vec := []float64{1888.62539593, 3303.82880589, 5179.86207441, 6355.6854363}
	runTestsZ0Z2(cos.ComovingDistanceZ1Z2Integrate, zLambdaCDM, exp_vec, distTol, t)
}

func TestLambdaCDMComovingDistanceZ1Z2Elliptic(t *testing.T) {
	cos := LambdaCDM{Om0: 0.3, Ol0: 0.7, H0: 70, Tcmb0: 0.}
	exp_vec := []float64{1888.62539593, 3303.82880589, 5179.86207441, 6355.6854363}
	runTestsZ0Z2(cos.ComovingDistanceZ1Z2Elliptic, zLambdaCDM, exp_vec, distTol, t)
}

func TestLambdaCDMLookbackTime(t *testing.T) {
	cos := LambdaCDM{Om0: 0.3, Ol0: 0.7, H0: 70, Tcmb0: 0.}
	// FlatLambdaCDM(70, 0.3).lookback_time
	exp_vec := []float64{5.04063793, 7.715337, 10.24035689, 11.35445676}
	runTests(cos.LookbackTime, zLambdaCDM, exp_vec, ageTol, t)
}

func TestLambdaCDMLookbackTimeIntegrate(t *testing.T) {
	cos := LambdaCDM{Om0: 0.3, Ol0: 0.7, H0: 70, Tcmb0: 0.}
	// FlatLambdaCDM(70, 0.3).lookback_time
	exp_vec := []float64{5.04063793, 7.715337, 10.24035689, 11.35445676}
	runTests(cos.LookbackTimeIntegrate, zLambdaCDM, exp_vec, ageTol, t)
}

func TestLambdaCDMLookbackTimeOM(t *testing.T) {
	cos := LambdaCDM{Om0: 0.3, Ol0: 0., H0: 70, Tcmb0: 0.}
	// FlatLambdaCDM(70, 0.3).lookback_time
	exp_vec := []float64{4.51471693, 6.62532254, 8.57486509, 9.45923582}
	runTests(cos.LookbackTimeOM, zLambdaCDM, exp_vec, ageTol, t)
}

func TestLambdaCDMLookbackTimeOL(t *testing.T) {
	cos := LambdaCDM{Om0: 0., Ol0: 0.5, H0: 70, Tcmb0: 0.}
	// FlatLambdaCDM(70, 0.3).lookback_time
	exp_vec := []float64{5.0616361, 7.90494991, 10.94241739, 12.52244605}
	runTests(cos.LookbackTimeOL, zLambdaCDM, exp_vec, ageTol, t)
}

func TestLambdaCDMAge(t *testing.T) {
	cos := LambdaCDM{Om0: 0.3, Ol0: 0.6, H0: 70, Tcmb0: 0.}
	// LambdaCDM(70, 0.3, 0.6).age(z)
	exp_vec := []float64{8.11137578, 5.54558439, 3.13456008, 2.06445301}

	runTests(cos.Age, zLambdaCDM, exp_vec, ageTol, t)
}

func TestLambdaCDMAgeFlatLCDM(t *testing.T) {
	cos := LambdaCDM{Om0: 0.3, Ol0: 0.7, H0: 70, Tcmb0: 0.}
	// FlatLambdaCDM(70, 0.3).age(z)
	exp_vec := []float64{8.42634602, 5.75164694, 3.22662706, 2.11252719}

	runTests(cos.Age, zLambdaCDM, exp_vec, ageTol, t)
}

func TestLambdaCDMAgeIntegrate(t *testing.T) {
	cos := LambdaCDM{Om0: 0.3, Ol0: 0.7, H0: 70, Tcmb0: 0.}
	// FlatLambdaCDM(70, 0.3).age(z)
	exp_vec := []float64{8.42634602, 5.75164694, 3.22662706, 2.11252719}

	runTests(cos.AgeIntegrate, zLambdaCDM, exp_vec, ageTol, t)
}

func TestLambdaCDMAgeOM(t *testing.T) {
	cos := LambdaCDM{Om0: 0.3, Ol0: 0., H0: 70, Tcmb0: 0.}
	// LambdaCDM(70, 0.3, 0.).age(z)
	exp_vec := []float64{6.78287955, 4.67227393, 2.72273139, 1.83836065}
	runTests(cos.AgeOM, zLambdaCDM, exp_vec, ageTol, t)
	runTests(cos.AgeIntegrate, zLambdaCDM, exp_vec, ageTol, t)
}

func TestLambdaCDMAgeOL(t *testing.T) {
	cos := LambdaCDM{Om0: 0., Ol0: 0.5, H0: 70, Tcmb0: 0.}
	// FlatLambdaCDM(70, 0.3).lookback_time
	exp_vec := []float64{12.34935796, 9.50604415, 6.46857667, 4.88854801}
	runTests(cos.AgeOL, zLambdaCDM, exp_vec, ageTol, t)
}

// Analytic case of Omega_Lambda = 0
func TestLambdaCDMEOm(t *testing.T) {
	zLambdaCDM := []float64{1.0, 10.0, 500.0, 1000.0}
	cos := LambdaCDM{Om0: 1.0, Ol0: 0., H0: 70, Tcmb0: 0.}
	hubbleDistance := SpeedOfLightKmS / cos.H0
	exp_vec := make([]float64, len(zLambdaCDM))
	for i, z := range zLambdaCDM {
		exp_vec[i] = 2.0 * hubbleDistance * (1 - math.Sqrt(1/(1+z)))
	}
	runTests(cos.ComovingDistance, zLambdaCDM, exp_vec, distTol, t)
}

func TestLambdaCDMEvecLcdm(t *testing.T) {
	cos := LambdaCDM{Om0: 0.27, Ol0: 0.73, H0: 70, Tcmb0: 0.}
	// Check array
	z := []float64{0.5, 1.0}
	// FlatLambdaCDM (v1.3.2)
	exp := []float64{1.2811127975318957, 1.7}
	obs := cos.Evec(z)
	if !floats.EqualApprox(obs, exp, eTol) {
		t.Errorf("Failed array float LCDM test.  Expected %v, return %v", exp, obs)
	}
}
