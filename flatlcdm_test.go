package cosmo

import (
	"math"
	"strings"
	"testing"
)

var zFlatLCDM = []float64{0.5, 1.0, 2.0, 3.0}

// Calculated via Python AstroPy
//   from astropy.cosmology import FlatLambdaCDM
//   z = np.asarray([0.5, 1.0, 2.0, 3.0])
var testTableFlatLCDM = map[string]struct {
	cos      FlatLCDM
	function string
	exp      []float64
}{
	"FlatLCDMDistanceModulus":    {FlatLCDM{H0: 70, Om0: 0.3}, "DistanceModulus", []float64{42.26118542, 44.10023766, 45.95719725, 47.02611193}},
	"FlatLCDMLuminosityDistance": {FlatLCDM{H0: 70, Om0: 0.3}, "LuminosityDistance", []float64{2832.9380939, 6607.65761177, 15539.58622323, 25422.74174519}},
	// Calculated via FlatLambdaCDM(70, 1.0).comoving_distance(z)
	"FlatLCDMComovingDistanceEdS":           {FlatLCDM{H0: 70, Om0: 1}, "ComovingDistance", []float64{1571.79831586, 2508.77651427, 3620.20576208, 4282.7494}},
	"FlatLCDMAngularDiameterDistance":       {FlatLCDM{H0: 70, Om0: 0.3}, "AngularDiameterDistance", []float64{1259.08359729, 1651.91440294, 1726.62069147, 1588.92135907}},
	"FlatLCDMComovingTransverseDistance":    {FlatLCDM{H0: 70, Om0: 0.3}, "ComovingTransverseDistance", []float64{1888.62539593, 3303.82880589, 5179.86207441, 6355.6854363}},
	"FlatLCDMComovingDistanceZ1Z2Integrate": {FlatLCDM{H0: 70, Om0: 0.3}, "ComovingDistanceZ1Z2", []float64{1888.62539593, 3303.82880589, 5179.86207441, 6355.6854363}},
	"FlatLCDMComovingDistanceZ1Z2Elliptic":  {FlatLCDM{H0: 70, Om0: 0.3}, "ComovingDistanceZ1Z2", []float64{1888.62539593, 3303.82880589, 5179.86207441, 6355.6854363}},
	// Calculated via FlatLambdaCDM(70, 0.3).lookback_time(z)
	"FlatLCDMLookbackTime": {FlatLCDM{H0: 70, Om0: 0.3}, "LookbackTime", []float64{5.04063793, 7.715337, 10.24035689, 11.35445676}},
	// Calculated via FlatLambdaCDM(70, 1.0).lookback_time(z)
	"FlatLCDMLookbackTimeEdS": {FlatLCDM{H0: 70, Om0: 1.0}, "LookbackTime", []float64{4.24332906, 6.0199092, 7.52015258, 8.14826851}},
	//   FlatLambdaCDM(70, 0.3).age(z)
	"FlatLCDMAge": {FlatLCDM{H0: 70, Om0: 0.3}, "Age", []float64{8.42634602, 5.75164694, 3.22662706, 2.11252719}},
	//   FlatLambdaCDM(70, 1.0).age(z)
	"FlatLCDMAgeEdS": {FlatLCDM{H0: 70, Om0: 1.0}, "Age", []float64{5.06897781, 3.29239767, 1.79215429, 1.16403836}},
}

func TestTableFlatLCDM(t *testing.T) {
	for _, test := range testTableFlatLCDM {
		switch {
		case strings.HasSuffix(test.function, "Z1Z2"):
			runTestsZ0Z2ByName(test.cos, test.function, zFlatLCDM, test.exp, distTol, t)
		default:
			runTestsByName(test.cos, test.function, zFlatLCDM, test.exp, distTol, t)
		}
	}
}

func TestFlatLCDMCosmologyInterface(t *testing.T) {
	ageDistance := func(cos FLRW) {
		z := 0.5
		age := cos.Age(z)
		dc := cos.ComovingDistance(z)
		_, _ = age, dc
	}

	cos := FlatLCDM{H0: 70, Om0: 0.27}
	ageDistance(cos)
}

// TestE* tests that basic calculation of E
//   https://github.com/astropy/astropy/blob/master/astropy/cosmology/tests/test_cosmology.py
func TestFlatLCDME(t *testing.T) {
	var z, exp float64
	cos := FlatLCDM{H0: 70, Om0: 0.27}

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

// Analytic case of Omega_Lambda = 0
func TestFlatLCDMEOm(t *testing.T) {
	cos := FlatLCDM{H0: 70, Om0: 1.0}
	zVec := []float64{1.0, 10.0, 500.0, 1000.0}
	hubbleDistance := SpeedOfLightKmS / cos.H0
	expVec := make([]float64, len(zVec))
	for i, z := range zVec {
		expVec[i] = 2.0 * hubbleDistance * (1 - math.Sqrt(1/(1+z)))
	}
	runTests(cos.ComovingDistance, zVec, expVec, distTol, t)
}
