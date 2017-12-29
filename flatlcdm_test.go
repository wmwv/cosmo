package cosmo

import (
	"math"
	"testing"
)

var zFlatLCDM = []float64{0.5, 1.0, 2.0, 3.0}

// Calculated via Python AstroPy
//   from astropy.cosmology import FlatLambdaCDM
//   z = np.asarray([0.5, 1.0, 2.0, 3.0])
var answersFlatLCDM = map[string][]float64{
	"FlatLCDMDistanceModulus":    []float64{42.26118542, 44.10023766, 45.95719725, 47.02611193},
	"FlatLCDMLuminosityDistance": []float64{2832.9380939, 6607.65761177, 15539.58622323, 25422.74174519},
	// Calculated via FlatLambdaCDM(70, 1.0).comoving_distance(z)
	"FlatLCDMComovingDistanceEdS":           []float64{1571.79831586, 2508.77651427, 3620.20576208, 4282.7494},
	"FlatLCDMAngularDiameterDistance":       []float64{1259.08359729, 1651.91440294, 1726.62069147, 1588.92135907},
	"FlatLCDMComovingTransverseDistance":    []float64{1888.62539593, 3303.82880589, 5179.86207441, 6355.6854363},
	"FlatLCDMComovingDistanceZ1Z2Integrate": []float64{1888.62539593, 3303.82880589, 5179.86207441, 6355.6854363},
	"FlatLCDMComovingDistanceZ1Z2Elliptic":  []float64{1888.62539593, 3303.82880589, 5179.86207441, 6355.6854363},
	// Calculated via FlatLambdaCDM(70, 0.3).lookback_time(z)
	"FlatLCDMLookbackTime": []float64{5.04063793, 7.715337, 10.24035689, 11.35445676},
	// Calculated via FlatLambdaCDM(70, 1.0).lookback_time(z)
	"FlatLCDMLookbackTimeEdS": []float64{4.24332906, 6.0199092, 7.52015258, 8.14826851},
	// Calculated via FlatLambdaCDM(70, 0.3).lookback_time(z)
	"FlatLCDMLookbackTimeIntegrate": []float64{5.04063793, 7.715337, 10.24035689, 11.35445676},
	//   FlatLambdaCDM(70, 0.3).age(z)
	"FlatLCDMAge": []float64{8.42634602, 5.75164694, 3.22662706, 2.11252719},
	//   FlatLambdaCDM(70, 1.0).age(z)
	"FlatLCDMAgeEdS": []float64{5.06897781, 3.29239767, 1.79215429, 1.16403836},
	//   FlatLambdaCDM(70, 0.3).age(z)
	"FlatLCDMAgeIntegrate": []float64{8.42634602, 5.75164694, 3.22662706, 2.11252719},
}

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
	var z, exp float64
	cos := FlatLCDM{Om0: 0.27, H0: 70, Tcmb0: 0.}

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

func TestFlatLCDMDistanceModulus(t *testing.T) {
	cos := FlatLCDM{Om0: 0.3, H0: 70, Tcmb0: 0.}
	exp_vec := answersFlatLCDM["FlatLCDMDistanceModulus"]
	runTests(cos.DistanceModulus, zFlatLCDM, exp_vec, distmodTol, t)
}

func TestFlatLCDMLuminosityDistance(t *testing.T) {
	cos := FlatLCDM{Om0: 0.3, H0: 70, Tcmb0: 0.}
	exp_vec := answersFlatLCDM["FlatLCDMLuminosityDistance"]
	runTests(cos.LuminosityDistance, zFlatLCDM, exp_vec, distmodTol, t)
}

func TestFlatLCDMAngularDiameterDistance(t *testing.T) {
	cos := FlatLCDM{Om0: 0.3, H0: 70, Tcmb0: 0.}
	exp_vec := answersFlatLCDM["FlatLCDMAngularDiameterDistance"]
	runTests(cos.AngularDiameterDistance, zFlatLCDM, exp_vec, distTol, t)
}

func TestFlatLCDMComovingTransverseDistance(t *testing.T) {
	cos := FlatLCDM{Om0: 0.3, H0: 70, Tcmb0: 0.}
	exp_vec := answersFlatLCDM["FlatLCDMComovingTransverseDistance"]
	runTests(cos.ComovingTransverseDistance, zFlatLCDM, exp_vec, distTol, t)
}

func TestFlatLCDMComovingDistanceZ1Z2Integrate(t *testing.T) {
	cos := FlatLCDM{Om0: 0.3, H0: 70, Tcmb0: 0.}
	exp_vec := answersFlatLCDM["FlatLCDMComovingDistanceZ1Z2Integrate"]
	runTestsZ0Z2(cos.comovingDistanceZ1Z2Integrate, zFlatLCDM, exp_vec, distTol, t)
}

func TestFlatLCDMComovingDistanceZ1Z2Elliptic(t *testing.T) {
	cos := FlatLCDM{Om0: 0.3, H0: 70, Tcmb0: 0.}
	exp_vec := answersFlatLCDM["FlatLCDMComovingDistanceZ1Z2Elliptic"]
	runTestsZ0Z2(cos.comovingDistanceZ1Z2Elliptic, zFlatLCDM, exp_vec, distTol, t)
}

func TestFlatLCDMLookbackTime(t *testing.T) {
	cos := FlatLCDM{Om0: 0.3, H0: 70, Tcmb0: 0.}
	exp_vec := answersFlatLCDM["FlatLCDMLookbackTime"]
	runTests(cos.LookbackTime, zFlatLCDM, exp_vec, ageTol, t)
}

func TestFlatLCDMLookbackTimeEdS(t *testing.T) {
	cos := FlatLCDM{Om0: 1.0, H0: 70, Tcmb0: 0.}
	exp_vec := answersFlatLCDM["FlatLCDMLookbackTimeEdS"]
	runTests(cos.LookbackTime, zFlatLCDM, exp_vec, ageTol, t)
}

func TestFlatLCDMLookbackTimeIntegrate(t *testing.T) {
	cos := FlatLCDM{Om0: 0.3, H0: 70, Tcmb0: 0.}
	exp_vec := answersFlatLCDM["FlatLCDMLookbackTimeIntegrate"]
	runTests(cos.lookbackTimeIntegrate, zFlatLCDM, exp_vec, ageTol, t)
}

func TestFlatLCDMAge(t *testing.T) {
	cos := FlatLCDM{Om0: 0.3, H0: 70, Tcmb0: 0.}
	exp_vec := answersFlatLCDM["FlatLCDMAge"]
	runTests(cos.Age, zFlatLCDM, exp_vec, ageTol, t)
	runTests(cos.ageIntegrate, zFlatLCDM, exp_vec, ageTol, t)
}

func TestFlatLCDMAgeIntegrate(t *testing.T) {
	cos := FlatLCDM{Om0: 0.3, H0: 70, Tcmb0: 0.}
	exp_vec := answersFlatLCDM["FlatLCDMAgeIntegrate"]
	runTests(cos.ageIntegrate, zFlatLCDM, exp_vec, ageTol, t)
}

// Analytic case of Omega_Lambda = 0
func TestFlatLCDMEOm(t *testing.T) {
	cos := FlatLCDM{Om0: 1.0, H0: 70, Tcmb0: 0.}
	z_vec := []float64{1.0, 10.0, 500.0, 1000.0}
	hubbleDistance := SpeedOfLightKmS / cos.H0
	exp_vec := make([]float64, len(z_vec))
	for i, z := range z_vec {
		exp_vec[i] = 2.0 * hubbleDistance * (1 - math.Sqrt(1/(1+z)))
	}
	runTests(cos.ComovingDistance, z_vec, exp_vec, distTol, t)
}
