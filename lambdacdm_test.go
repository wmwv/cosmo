package cosmo

import (
	"math"
	"testing"
)

var zLambdaCDM = []float64{0.5, 1.0, 2.0, 3.0}

// Calculated via Python AstroPy
//   from astropy.cosmology import LambdaCDM
//   z = np.asarray([0.5, 1.0, 2.0, 3.0])
var testTableLambdaCDM = map[string]struct {
	cos      LambdaCDM
	function string
	exp      []float64
}{
	"LambdaCDMDistanceModulus": {LambdaCDM{H0: 70, Om0: 0.3, Ol0: 0.7}, "DistanceModulus", []float64{42.26118542, 44.10023766, 45.95719725, 47.02611193}},
	//   LambdaCDM(70, 0.3, 0.7).luminosity_distance(z)
	"LambdaCDMLuminosityDistanceFlat": {LambdaCDM{H0: 70, Om0: 0.3, Ol0: 0.7}, "LuminosityDistance", []float64{2832.9380939, 6607.65761177, 15539.58622323, 25422.74174519}},
	//   LambdaCDM(70, 0.3, 0.6).luminosity_distance(z)
	"LambdaCDMLuminosityDistanceNonflat":     {LambdaCDM{H0: 70, Om0: 0.3, Ol0: 0.6}, "LuminosityDistance", []float64{2787.51504671, 6479.83450953, 15347.21516211, 25369.7240234}},
	"LambdaCDMAngularDiameterDistance":       {LambdaCDM{H0: 70, Om0: 0.3, Ol0: 0.7}, "AngularDiameterDistance", []float64{1259.08359729, 1651.91440294, 1726.62069147, 1588.92135907}},
	"LambdaCDMComovingTransverseDistance":    {LambdaCDM{H0: 70, Om0: 0.3, Ol0: 0.7}, "ComovingTransverseDistance", []float64{1888.62539593, 3303.82880589, 5179.86207441, 6355.6854363}},
	"LambdaCDMComovingDistanceZ1Z2Integrate": {LambdaCDM{H0: 70, Om0: 0.3, Ol0: 0.7}, "ComovingDistanceZ1Z2Integrate", []float64{1888.62539593, 3303.82880589, 5179.86207441, 6355.6854363}},
	"LambdaCDMComovingDistanceZ1Z2Elliptic":  {LambdaCDM{H0: 70, Om0: 0.3, Ol0: 0.7}, "ComovingDistanceZ1Z2Elliptic", []float64{1888.62539593, 3303.82880589, 5179.86207441, 6355.6854363}},
	// LambdaCDM(70, 0.3, 0).comoving_distance(z)
	"LambdaCDMComovingDistanceNonflatOM": {LambdaCDM{H0: 70, Om0: 0.3, Ol0: 0.}, "ComovingDistance", []float64{1679.81156606, 2795.15602075, 4244.25192263, 5178.38877021}},
	// LambdaCDM(70, 0.3, 0).comoving_transverse_distance(z)
	"LambdaCDMComovingTransverseDistanceNonflatOM": {LambdaCDM{H0: 70, Om0: 0.3, Ol0: 0.}, "ComovingTransverseDistance", []float64{1710.1240353, 2936.1472205, 4747.54480615, 6107.95517311}},
	// FlatLambdaCDM(70, 1.0).comoving_distance(z)
	"LambdaCDMComovingDistanceEdS": {LambdaCDM{H0: 70, Om0: 1.0, Ol0: 0.}, "ComovingDistance", []float64{1571.79831586, 2508.77651427, 3620.20576208, 4282.7494}},
	// LambdaCDM(70, 0.3, 0).lookback_time(z)
	"LambdaCDMLookbackTime": {LambdaCDM{H0: 70, Om0: 0.3, Ol0: 0.7}, "LookbackTime", []float64{5.04063793, 7.715337, 10.24035689, 11.35445676}},
	// LambdaCDM(70, 0.3, 0).lookback_time(z)
	"LambdaCDMLookbackTimeOM": {LambdaCDM{H0: 70, Om0: 0.3, Ol0: 0.}, "LookbackTime", []float64{4.51471693, 6.62532254, 8.57486509, 9.45923582}},
	// LambdaCDM(70, 0.3, 0.7).lookback_time(z)
	"LambdaCDMLookbackTimeOL": {LambdaCDM{H0: 70, Om0: 0., Ol0: 0.5}, "LookbackTime", []float64{5.0616361, 7.90494991, 10.94241739, 12.52244605}},
	"LambdaCDMAge":            {LambdaCDM{H0: 70, Om0: 0.3, Ol0: 0.6}, "Age", []float64{8.11137578, 5.54558439, 3.13456008, 2.06445301}},
	"LambdaCDMAgeFlatLCDM":    {LambdaCDM{H0: 70, Om0: 0.3, Ol0: 0.7}, "Age", []float64{8.42634602, 5.75164694, 3.22662706, 2.11252719}},
	// FlatLambdaCDM(70, 1.0).age(z)
	"LambdaCDMAgeEdS": {LambdaCDM{H0: 70, Om0: 1.0, Ol0: 0.}, "Age", []float64{5.06897781, 3.29239767, 1.79215429, 1.16403836}},
	// LambdaCDM(70, 0.3, 0.).age(z)
	"LambdaCDMAgeOM": {LambdaCDM{H0: 70, Om0: 0.3, Ol0: 0.}, "Age", []float64{6.78287955, 4.67227393, 2.72273139, 1.83836065}},
	// FlatLambdaCDM(70, 0, 0.5).lookback_time
	"LambdaCDMAgeOL": {LambdaCDM{H0: 70, Om0: 0., Ol0: 0.5}, "Age", []float64{12.34935796, 9.50604415, 6.46857667, 4.88854801}},
}

func TestLambdaCDMCosmologyInterface(t *testing.T) {
	age_distance := func(cos FLRW) {
		z := 0.5
		age := cos.Age(z)
		dc := cos.ComovingDistance(z)
		_, _ = age, dc
	}

	cos := LambdaCDM{H0: 70, Om0: 0.27, Ol0: 0.73}
	age_distance(cos)
}

// TestE* tests that basic calculation of E
//   https://github.com/astropy/astropy/blob/master/astropy/cosmology/tests/test_cosmology.py
func TestLambdaCDMELcdm(t *testing.T) {
	var z, exp, tol float64
	cos := LambdaCDM{H0: 70, Om0: 0.27, Ol0: 0.73}

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
	test := testTableLambdaCDM["LambdaCDMDistanceModulus"]
	runTestsByName(test.cos, test.function, zLambdaCDM, test.exp, distTol, t)
}

func TestLambdaCDMLuminosityDistanceFlat(t *testing.T) {
	test := testTableLambdaCDM["LambdaCDMLuminosityDistanceFlat"]
	runTestsByName(test.cos, test.function, zLambdaCDM, test.exp, distTol, t)
}

func TestLambdaCDMLuminosityDistanceNonflat(t *testing.T) {
	test := testTableLambdaCDM["LambdaCDMLuminosityDistanceNonflat"]
	runTestsByName(test.cos, test.function, zLambdaCDM, test.exp, distTol, t)
}

func TestLambdaCDMAngularDiameterDistance(t *testing.T) {
	test := testTableLambdaCDM["LambdaCDMAngularDiameterDistance"]
	runTestsByName(test.cos, test.function, zLambdaCDM, test.exp, distTol, t)
}

func TestLambdaCDMComovingTransverseDistance(t *testing.T) {
	test := testTableLambdaCDM["LambdaCDMComovingTransverseDistance"]
	runTestsByName(test.cos, test.function, zLambdaCDM, test.exp, distTol, t)
}

func TestLambdaCDMComovingDistanceZ1Z2Integrate(t *testing.T) {
	test := testTableLambdaCDM["LambdaCDMComovingDistanceZ1Z2Integrate"]
	runTestsZ0Z2(test.cos.comovingDistanceZ1Z2Integrate, zLambdaCDM, test.exp, distTol, t)
}

func TestLambdaCDMComovingDistanceZ1Z2Elliptic(t *testing.T) {
	test := testTableLambdaCDM["LambdaCDMComovingDistanceZ1Z2Elliptic"]
	runTestsZ0Z2(test.cos.comovingDistanceZ1Z2Elliptic, zLambdaCDM, test.exp, distTol, t)
}

func TestLambdaCDMComovingDistanceNonflatOM(t *testing.T) {
	test := testTableLambdaCDM["LambdaCDMComovingDistanceNonflatOM"]
	runTestsByName(test.cos, test.function, zLambdaCDM, test.exp, distTol, t)
}

func TestLambdaCDMComovingTransverseDistanceNonflatOM(t *testing.T) {
	test := testTableLambdaCDM["LambdaCDMComovingTransverseDistanceNonflatOM"]
	runTestsByName(test.cos, test.function, zLambdaCDM, test.exp, distTol, t)
}

func TestLambdaCDMComovingDistanceEdS(t *testing.T) {
	test := testTableLambdaCDM["LambdaCDMComovingDistanceEdS"]
	runTestsByName(test.cos, test.function, zLambdaCDM, test.exp, distTol, t)
}

func TestLambdaCDMLookbackTime(t *testing.T) {
	test := testTableLambdaCDM["LambdaCDMLookbackTime"]
	runTestsByName(test.cos, test.function, zLambdaCDM, test.exp, ageTol, t)
}

func TestLambdaCDMLookbackTimeOM(t *testing.T) {
	test := testTableLambdaCDM["LambdaCDMLookbackTimeOM"]
	runTestsByName(test.cos, test.function, zLambdaCDM, test.exp, ageTol, t)
}

func TestLambdaCDMLookbackTimeOL(t *testing.T) {
	test := testTableLambdaCDM["LambdaCDMLookbackTimeOL"]
	runTestsByName(test.cos, test.function, zLambdaCDM, test.exp, ageTol, t)
}

func TestLambdaCDMAge(t *testing.T) {
	// LambdaCDM(70, 0.3, 0.6).age(z)
	test := testTableLambdaCDM["LambdaCDMAge"]
	runTestsByName(test.cos, test.function, zLambdaCDM, test.exp, ageTol, t)
}

func TestLambdaCDMAgeFlatLCDM(t *testing.T) {
	// FlatLambdaCDM(70, 0.3).age(z)
	test := testTableLambdaCDM["LambdaCDMAgeFlatLCDM"]
	runTestsByName(test.cos, test.function, zLambdaCDM, test.exp, ageTol, t)
}

func TestLambdaCDMAgeOM(t *testing.T) {
	test := testTableLambdaCDM["LambdaCDMAgeOM"]
	runTestsByName(test.cos, test.function, zLambdaCDM, test.exp, ageTol, t)
}

func TestLambdaCDMAgeEdS(t *testing.T) {
	test := testTableLambdaCDM["LambdaCDMAgeEdS"]
	runTestsByName(test.cos, test.function, zLambdaCDM, test.exp, ageTol, t)
}

func TestLambdaCDMAgeOL(t *testing.T) {
	test := testTableLambdaCDM["LambdaCDMAgeOL"]
	runTestsByName(test.cos, test.function, zLambdaCDM, test.exp, ageTol, t)
}

// Analytic case of Omega_Lambda = 0
func TestLambdaCDMEOm(t *testing.T) {
	zLambdaCDM := []float64{1.0, 10.0, 500.0, 1000.0}
	cos := LambdaCDM{H0: 70, Om0: 1.0, Ol0: 0.}
	hubbleDistance := SpeedOfLightKmS / cos.H0
	exp_vec := make([]float64, len(zLambdaCDM))
	for i, z := range zLambdaCDM {
		exp_vec[i] = 2.0 * hubbleDistance * (1 - math.Sqrt(1/(1+z)))
	}
	runTests(cos.ComovingDistance, zLambdaCDM, exp_vec, distTol, t)
}
