package cosmo

import (
	"math"
	"strings"
	"testing"
)

var zWCDM = []float64{0.5, 1.0, 2.0, 3.0}

// Calculated via Python AstroPy
//   from astropy.cosmology import wCDM
//   z = np.asarray([0.5, 1.0, 2.0, 3.0])
var testTableWCDM = map[string]struct {
	cos      WCDM
	function string
	exp      []float64
}{
	//   wCDM(70, 0.3, 0.7, -  "WCDMELcdm": WCDM{H0: 70, Om0: 0.27, Ol0: 0.73, W0: -1},
	"WCDMDistanceModulus":            {WCDM{H0: 70, Om0: 0.3, Ol0: 0.7, W0: -1.2}, "DistanceModulus", []float64{42.32710911, 44.17957201, 46.03118144, 47.09228735}},
	"WCDMLuminosityDistanceFlatLCDM": {WCDM{H0: 70, Om0: 0.3, Ol0: 0.7, W0: -1}, "LuminosityDistance", []float64{2832.9380939, 6607.65761177, 15539.58622323, 25422.74174519}},
	//   wCDM(70, 0.3, 0.9, -1).luminosity_distance(z)
	"WCDMLuminosityDistanceNegativeOkLCDM": {WCDM{H0: 70, Om0: 0.3, Ol0: 0.9, W0: -1}, "LuminosityDistance", []float64{2933.96568944, 6896.93040403, 15899.60122012, 25287.53295915}},
	//   wCDM(70, 0.3, 0.7, -1.1).luminosity_distance(z)
	"WCDMLuminosityDistanceFlat": {WCDM{H0: 70, Om0: 0.3, Ol0: 0.7, W0: -1.1}, "LuminosityDistance", []float64{2877.10314183, 6734.38177991, 15823.59621899, 25841.56448508}},
	//   wCDM(70, 0.3, 0.6, -0.8).luminosity_distance(z)
	"WCDMLuminosityDistanceNonflat":     {WCDM{H0: 70, Om0: 0.3, Ol0: 0.6, W0: -0.8}, "LuminosityDistance", []float64{2713.4660301, 6257.24866642, 14794.59911147, 24496.30592953}},
	"WCDMAngularDiameterDistance":       {WCDM{H0: 70, Om0: 0.3, Ol0: 0.7, W0: -1}, "AngularDiameterDistance", []float64{1259.08359729, 1651.91440294, 1726.62069147, 1588.92135907}},
	"WCDMComovingTransverseDistance":    {WCDM{H0: 70, Om0: 0.3, Ol0: 0.7, W0: -1}, "ComovingTransverseDistance", []float64{1888.62539593, 3303.82880589, 5179.86207441, 6355.6854363}},
	"WCDMComovingDistanceZ1Z2Integrate": {WCDM{H0: 70, Om0: 0.3, Ol0: 0.7, W0: -1}, "ComovingDistanceZ1Z2", []float64{1888.62539593, 3303.82880589, 5179.86207441, 6355.6854363}},
	"WCDMComovingDistanceZ1Z2Elliptic":  {WCDM{H0: 70, Om0: 0.3, Ol0: 0.7, W0: -1}, "ComovingDistanceZ1Z2", []float64{1888.62539593, 3303.82880589, 5179.86207441, 6355.6854363}},
	//   wCDM(70, 1.0, 0., -1, 0).comoving_distance(z)
	"WCDMComovingDistanceEdS": {WCDM{H0: 70, Om0: 1.0, Ol0: 0., W0: -1}, "ComovingDistance", []float64{1571.79831586, 2508.77651427, 3620.20576208, 4282.7494}},
	//   wCDM(70, 0.3, 0.0, -1).comoving_distance(z)
	"WCDMComovingDistanceNonflatOM": {WCDM{H0: 70, Om0: 0.3, Ol0: 0., W0: -1}, "ComovingDistance", []float64{1679.81156606, 2795.15602075, 4244.25192263, 5178.38877021}},
	//   wCDM(70, 0.3, 0.0, -1).comoving_transverse_distance(z)
	"WCDMComovingTransverseDistanceNonflatOM": {WCDM{H0: 70, Om0: 0.3, Ol0: 0., W0: -1}, "ComovingTransverseDistance", []float64{1710.1240353, 2936.1472205, 4747.54480615, 6107.95517311}},
	// Calculated via astropy.cosmology.wCDM(70, 0.3, 0.7, w0=-1.2).lookback_time
	"WCDMLookbackTime": {WCDM{H0: 70, Om0: 0.3, Ol0: 0.7, W0: -1.2}, "LookbackTime", []float64{5.18796426, 7.98542226, 10.58842012, 11.71902479}},
	// Calculated via astropy.cosmology.wCDM(70, 0.3, 0.7, w0=-1.1).lookback_time
	"WCDMLookbackTimeIntegrate": {WCDM{H0: 70, Om0: 0.3, Ol0: 0.7, W0: -1.1}, "LookbackTime", []float64{5.11509518, 7.85406053, 10.42213038, 11.54588106}},
	// Calculated via astropy.cosmology.wCDM(70, 0.3, 0.7, w0=-0.9).lookback_time
	"WCDMLookbackTimeOM": {WCDM{H0: 70, Om0: 0.3, Ol0: 0., W0: -0.9}, "LookbackTime", []float64{4.51471693, 6.62532254, 8.57486509, 9.45923582}},
	// Calculated via astropy.cosmology.wCDM(70, 0, 0.5, -0.9).lookback_time
	"WCDMLookbackTimeOL": {WCDM{H0: 70, Om0: 0, Ol0: 0.5, W0: -0.9}, "LookbackTime", []float64{5.00576631, 7.78841245, 10.76147941, 12.31462586}},
	//   astropy.cosmology.WCDM(70, 0.3, 0.6).age(z)
	"WCDMAge": {WCDM{H0: 70, Om0: 0.3, Ol0: 0.6, W0: -1}, "Age", []float64{8.11137578, 5.54558439, 3.13456008, 2.06445301}},
	//   astropy.cosmology.wCDM(70, 0.3).age(z)
	"WCDMAgeFlatLCDM": {WCDM{H0: 70, Om0: 0.3, Ol0: 0.7, W0: -1}, "Age", []float64{8.42634602, 5.75164694, 3.22662706, 2.11252719}},
	//   astropy.cosmology.wCDM(70, 0.3).age(z)
	"WCDMAgeIntegrate": {WCDM{H0: 70, Om0: 0.3, Ol0: 0.7, W0: -1}, "Age", []float64{8.42634602, 5.75164694, 3.22662706, 2.11252719}},
	//   astropy.cosmology.WCDM(70, 0.3, 0.).age(z)
	"WCDMAgeOM": {WCDM{H0: 70, Om0: 0.3, Ol0: 0., W0: -1}, "Age", []float64{6.78287955, 4.67227393, 2.72273139, 1.83836065}},
	//   LambdaCDM(70, 0, 0.5).age(z)
	"WCDMAgeOL": {WCDM{H0: 70, Om0: 0, Ol0: 0.5, W0: -1}, "Age", []float64{12.34935796, 9.50604415, 6.46857667, 4.88854801}},
}

func TestTableWCDM(t *testing.T) {
	for _, test := range testTableWCDM {
		switch {
		case strings.HasSuffix(test.function, "Z1Z2"):
			runTestsZ0Z2ByName(test.cos, test.function, zWCDM, test.exp, distTol, t)
		default:
			runTestsByName(test.cos, test.function, zWCDM, test.exp, distTol, t)
		}
	}
}

func TestWCDMCosmologyInterface(t *testing.T) {
	age_distance := func(cos FLRW) {
		z := 0.5
		age := cos.Age(z)
		dc := cos.ComovingDistance(z)
		_, _ = age, dc
	}

	cos := WCDM{H0: 70, Om0: 0.27, Ol0: 0.73, W0: -1}
	age_distance(cos)
}

// TestE* tests that basic calculation of E
//   https://github.com/astropy/astropy/blob/master/astropy/cosmology/tests/test_cosmology.py
func TestWCDMELcdm(t *testing.T) {
	cos := WCDM{H0: 70, Om0: 0.27, Ol0: 0.73, W0: -1}
	// Check value of E(z=1.0)
	//   OM, OL, OK, z = 0.27, 0.73, 0.0, 1.0
	//   sqrt(OM*(1+z)**3 + OK * (1+z)**2 + OL)
	//   sqrt(0.27*(1+1.0)**3 + 0.0 * (1+1.0)**2 + 0.73)
	//   sqrt(0.27*8 + 0 + 0.73)
	//   sqrt(2.89)
	z := 1.0
	exp := 1.7
	runTest(cos.E, z, exp, eTol, t, 0)

	exp = 1 / 1.7
	runTest(cos.Einv, z, exp, eTol, t, 0)
}

// Analytic case of Omega_Lambda = 0
func TestWCDMEOm(t *testing.T) {
	cos := WCDM{H0: 70, Om0: 1.0, Ol0: 0., W0: -1}
	z_vec := []float64{1.0, 10.0, 500.0, 1000.0}
	exp_vec := make([]float64, len(z_vec))
	hubbleDistance := SpeedOfLightKmS / cos.H0
	for i, z := range z_vec {
		exp_vec[i] = 2.0 * hubbleDistance * (1 - math.Sqrt(1/(1+z)))
	}
	runTests(cos.ComovingDistance, z_vec, exp_vec, distTol, t)
}
