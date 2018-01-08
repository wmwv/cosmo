package cosmo

import (
	"strings"
	"testing"
)

var zWACDM = []float64{0.5, 1.0, 2.0, 3.0}

// Calculated via Python AstroPy
//   from astropy.cosmology import w0waCDM
//   z = np.asarray([0.5, 1.0, 2.0, 3.0])
var testTableWACDM = map[string]struct {
	cos      WACDM
	function string
	exp      []float64
}{
	//   w0waCDM(70, 0.3, 0.7, -1.2, 3).distmod(z)
	"WACDMDistanceModulus": {WACDM{H0: 70, Om0: 0.3, Ol0: 0.7, W0: -1.2, WA: 3}, "DistanceModulus", []float64{42.20567831, 43.92122272, 45.57180818, 46.47483095}},
	//   w0waCDM(70, 0.3, 0.7, -1.2, 0).luminosity_distance(z)
	"WACDMLcdmDistanceModulus": {WACDM{H0: 70, Om0: 0.3, Ol0: 0.7, W0: -1.2, WA: 0}, "DistanceModulus", []float64{42.32710911, 44.17957201, 46.03118144, 47.09228735}},
	//   w0waCDM(70, 0.3, 0.7, -0.9, 2).luminosity_distance(z)
	"WACDMLuminosityDistanceFlat": {WACDM{H0: 70, Om0: 0.3, Ol0: 0.7, W0: -0.9, WA: 2}, "LuminosityDistance", []float64{2676.62203931, 5904.08905744, 12867.17142278, 19961.9490794}},
	//   w0waCDM(70, 0.3, 0.6, -1, 0).luminosity_distance(z)
	"WACDMLuminosityDistancePositiveOkLCDM": {WACDM{H0: 70, Om0: 0.3, Ol0: 0.6, W0: -1, WA: 0}, "LuminosityDistance", []float64{2787.51504671, 6479.83450953, 15347.21516211, 25369.7240234}},
	//   w0waCDM(70, 0.3, 0.9, -1, 0).luminosity_distance(z)
	"WACDMLuminosityDistanceNegativeOkLCDM": {WACDM{H0: 70, Om0: 0.3, Ol0: 0.9, W0: -1, WA: 0}, "LuminosityDistance", []float64{2933.96568944, 6896.93040403, 15899.60122012, 25287.53295915}},
	//   w0waCDM(70, 0.3, 0.6, -0.9, 2).luminosity_distance(z)
	"WACDMLuminosityDistanceNonflat": {WACDM{H0: 70, Om0: 0.3, Ol0: 0.6, W0: -0.9, WA: 2}, "LuminosityDistance", []float64{2659.67537448, 5901.12663329, 13049.93089016, 20468.18548013}},
	//   w0waCDM(70, 0.3, 0.7, -0.8, 2.5).angular_diameter_distance(z)
	"WACDMAngularDiameterDistance": {WACDM{H0: 70, Om0: 0.3, Ol0: 0.7, W0: -0.8, WA: 2.5}, "AngularDiameterDistance", []float64{1155.52181127, 1393.61319898, 1282.08090454, 1073.63096224}},
	//   w0waCDM(70, 0.3, 0., -1, 0).comoving_distance(z)
	"WACDMComovingDistanceNonflatOM": {WACDM{H0: 70, Om0: 0.3, Ol0: 0, W0: -1, WA: 0}, "ComovingDistance", []float64{1679.81156606, 2795.15602075, 4244.25192263, 5178.38877021}},
	//   w0waCDM(70, 0.3, 0., -1, 0).comoving_transverse_distance(z)
	"WACDMComovingTransverseDistanceNonflatOM": {WACDM{H0: 70, Om0: 0.3, Ol0: 0, W0: -1, WA: 0}, "ComovingTransverseDistance", []float64{1710.1240353, 2936.1472205, 4747.54480615, 6107.95517311}},
	//   w0waCDM(70, 1.0, 0., -1, 0).comoving_transverse_distance(z)
	"WACDMComovingDistanceEdS": {WACDM{H0: 70, Om0: 1.0, Ol0: 0, W0: -1, WA: 0}, "ComovingDistance", []float64{1571.79831586, 2508.77651427, 3620.20576208, 4282.7494}},
	//   w0waCDM(70, 0.3, 0.7, -1.2, -1.2).comoving_transverse_distance(z)
	"WACDMComovingTransverseDistance": {WACDM{H0: 70, Om0: 0.3, Ol0: 0.7, W0: -1.2, WA: -1.2}, "ComovingTransverseDistance", []float64{1985.54631561, 3533.91345688, 5524.66720808, 6731.56420461}},
	//   w0waCDM(70, 0.3, 0.6, -1, 0).comoving_transverse_distance(z)
	"WACDMComovingTransverseDistancePositiveOkLCDM": {WACDM{H0: 70, Om0: 0.3, Ol0: 0.6, W0: -1, WA: 0}, "ComovingTransverseDistance", []float64{1858.34336447, 3239.91725476, 5115.73838737, 6342.43100585}},
	//   w0waCDM(70, 0.3, 0.9, -1, 0).comoving_transverse_distance(z)
	"WACDMComovingTransverseDistanceNegativeOkLCDM": {WACDM{H0: 70, Om0: 0.3, Ol0: 0.9, W0: -1, WA: 0}, "ComovingTransverseDistance", []float64{1955.97712629, 3448.46520202, 5299.86707337, 6321.88323979}},
	//   w0waCDM(70, 0.3, 0.7, -0.9, 3.5)._comoving_distance_z1z2(0, z)
	"WACDMComovingDistanceZ1Z2Integrate": {WACDM{H0: 70, Om0: 0.3, Ol0: 0.7, W0: -0.9, WA: 3.5}, "ComovingDistanceZ1Z2", []float64{1726.71519955, 2709.17698433, 3538.63486291, 3798.28908226}},
	//   w0waCDM(70, 0.3, 0.7, -0.9, 3.5).lookback_time(z)
	"WACDMLookbackTime": {WACDM{H0: 70, Om0: 0.3, Ol0: 0.7, W0: -0.9, WA: 3.5}, "LookbackTime", []float64{4.64427098, 6.51439755, 7.65559243, 7.90553458}},
	//   w0waCDM(70, 0.3, 0.7, -1.1, 2.8).lookback_time(z)
	"WACDMLookbackTimeIntegrate": {WACDM{H0: 70, Om0: 0.3, Ol0: 0.7, W0: -1.1, WA: 2.8}, "LookbackTime", []float64{4.86957219, 7.09750717, 8.81344126, 9.36817438}},
	//   w0waCDM(70, 0.3, 0., -0.9, 2).lookback_time(z)
	"WACDMLookbackTimeOM": {WACDM{H0: 70, Om0: 0.3, Ol0: 0, W0: -0.9, WA: 2}, "LookbackTime", []float64{4.51471693, 6.62532254, 8.57486509, 9.45923582}},
	//   w0waCDM(70, 0., 0.5, -1, 0).lookback_time(z)
	"WACDMLookbackTimeOL": {WACDM{H0: 70, Om0: 0, Ol0: 0.5, W0: -1, WA: 0}, "LookbackTime", []float64{5.0616361, 7.90494991, 10.94241739, 12.52244605}},
	//   w0waCDM(70, 0.3, 0.6, -0.6, 3.5).age(z)
	"WACDMAge": {WACDM{H0: 70, Om0: 0.3, Ol0: 0.6, W0: -0.6, WA: 3.5}, "Age", []float64{2.70980463, 1.08619498, 0.21688951, 0.058307}},
	//   LambdaCDM(70, 0.3, 0.).age(z)
	"WACDMAgeOM": {WACDM{H0: 70, Om0: 0.3, Ol0: 0, W0: -1, WA: 0}, "Age", []float64{6.78287955, 4.67227393, 2.72273139, 1.83836065}},
	//   LambdaCDM(70, 0, 0.5).age(z)
	"WACDMAgeOL": {WACDM{H0: 70, Om0: 0, Ol0: 0.5, W0: -1, WA: 0}, "Age", []float64{12.34935796, 9.50604415, 6.46857667, 4.88854801}},
}

func TestTableWACDM(t *testing.T) {
	for _, test := range testTableWACDM {
		switch {
		case strings.HasSuffix(test.function, "Z1Z2"):
			runTestsZ0Z2ByName(test.cos, test.function, zWACDM, test.exp, distTol, t)
		default:
			runTestsByName(test.cos, test.function, zWACDM, test.exp, distTol, t)
		}
	}
}

func TestWACDMCosmologyInterface(t *testing.T) {
	age_distance := func(cos FLRW) {
		z := 0.5
		age := cos.Age(z)
		dc := cos.ComovingDistance(z)
		_, _ = age, dc
	}

	cos := WACDM{H0: 70, Om0: 0.27, Ol0: 0.73, W0: -1}
	age_distance(cos)
}

// TestE* tests that basic calculation of E
//   https://github.com/astropy/astropy/blob/master/astropy/cosmology/tests/test_cosmology.py
func TestWACDMELcdm(t *testing.T) {
	cos := WACDM{H0: 70, Om0: 0.27, Ol0: 0.73, W0: -1}
	var z, exp float64

	// Check value of E(z=1.0)
	//   OM, OL, OK, z = 0.27, 0.73, 0.0, 1.0
	//   sqrt(OM*(1+z)**3 + OK * (1+z)**2 + OL)
	//   sqrt(0.27*(1+1.0)**3 + 0.0 * (1+1.0)**2 + 0.73)
	//   sqrt(0.27*8 + 0 + 0.73)
	//   sqrt(2.89)
	z = 1.0
	exp = 1.7
	runTest(cos.E, z, exp, eTol, t, 0)

	exp = 1 / 1.7
	runTest(cos.Einv, z, exp, eTol, t, 0)
}
