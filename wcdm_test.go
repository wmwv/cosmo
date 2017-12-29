package cosmo

import (
	"math"
	"testing"
)

var zWCDM = []float64{0.5, 1.0, 2.0, 3.0}

// Calculated via Python AstroPy
//   from astropy.cosmology import w0waCDM
//   z = np.asarray([0.5, 1.0, 2.0, 3.0])
var answersWCDM = map[string][]float64{
	//   wCDM(70, 0.3, 0.7, -1.2).distmod(z)
	"WCDMDistanceModulus":           []float64{42.32710911, 44.17957201, 46.03118144, 47.09228735},
	"WCDMLuminosityDistanceFlatCDM": []float64{2832.9380939, 6607.65761177, 15539.58622323, 25422.74174519},
	//   wCDM(70, 0.3, 0.7, -1.1).luminosity_distance(z)
	"WCDMLuminosityDistanceFlat": []float64{2877.10314183, 6734.38177991, 15823.59621899, 25841.56448508},
	//   wCDM(70, 0.3, 0.6, -0.8).luminosity_distance(z)
	"WCDMLuminosityDistanceNonflat":     []float64{2713.4660301, 6257.24866642, 14794.59911147, 24496.30592953},
	"WCDMAngularDiameterDistance":       []float64{1259.08359729, 1651.91440294, 1726.62069147, 1588.92135907},
	"WCDMComovingTransverseDistance":    []float64{1888.62539593, 3303.82880589, 5179.86207441, 6355.6854363},
	"WCDMComovingDistanceZ1Z2Integrate": []float64{1888.62539593, 3303.82880589, 5179.86207441, 6355.6854363},
	"WCDMComovingDistanceZ1Z2Elliptic":  []float64{1888.62539593, 3303.82880589, 5179.86207441, 6355.6854363},
	//   wCDM(70, 1.0, 0., -1, 0).comoving_distance(z)
	"WCDMComovingDistanceEdS": []float64{1571.79831586, 2508.77651427, 3620.20576208, 4282.7494},
	//   wCDM(70, 0.3, 0.0, -1).comoving_distance(z)
	"WCDMComovingDistanceNonflatOM": []float64{1679.81156606, 2795.15602075, 4244.25192263, 5178.38877021},
	//   wCDM(70, 0.3, 0.0, -1).comoving_transverse_distance(z)
	"WCDMComovingTransverseDistanceNonflatOM": []float64{1710.1240353, 2936.1472205, 4747.54480615, 6107.95517311},
	// Calculated via astropy.cosmology.wCDM(70, 0.3, 0.7, w0=-1.2).lookback_time
	"WCDMLookbackTime": []float64{5.18796426, 7.98542226, 10.58842012, 11.71902479},
	// Calculated via astropy.cosmology.wCDM(70, 0.3, 0.7, w0=-1.1).lookback_time
	"WCDMLookbackTimeIntegrate": []float64{5.11509518, 7.85406053, 10.42213038, 11.54588106},
	// Calculated via astropy.cosmology.wCDM(70, 0.3, 0.7, w0=-0.9).lookback_time
	"WCDMLookbackTimeOM": []float64{4.51471693, 6.62532254, 8.57486509, 9.45923582},
	// Calculated via astropy.cosmology.wCDM(70, 0.3).lookback_time
	"WCDMLookbackTimeOL": []float64{5.0616361, 7.90494991, 10.94241739, 12.52244605},
	//   astropy.cosmology.WCDM(70, 0.3, 0.6).age(z)
	"WCDMAge": []float64{8.11137578, 5.54558439, 3.13456008, 2.06445301},
	//   astropy.cosmology.wCDM(70, 0.3).age(z)
	"WCDMAgeFlatLCDM": []float64{8.42634602, 5.75164694, 3.22662706, 2.11252719},
	//   astropy.cosmology.wCDM(70, 0.3).age(z)
	"WCDMAgeIntegrate": []float64{8.42634602, 5.75164694, 3.22662706, 2.11252719},
	//   astropy.cosmology.WCDM(70, 0.3, 0.).age(z)
	"WCDMAgeOM": []float64{6.78287955, 4.67227393, 2.72273139, 1.83836065},
	//   LambdaCDM(70, 0, 0.5).age(z)
	"WCDMAgeOL": []float64{12.34935796, 9.50604415, 6.46857667, 4.88854801},
}

func TestWCDMCosmologyInterface(t *testing.T) {
	age_distance := func(cos FLRW) {
		z := 0.5
		age := cos.Age(z)
		dc := cos.ComovingDistance(z)
		_, _ = age, dc
	}

	cos := WCDM{Om0: 0.27, Ol0: 0.73, W0: -1, H0: 70, Tcmb0: 0.}
	age_distance(cos)
}

// TestE* tests that basic calculation of E
//   https://github.com/astropy/astropy/blob/master/astropy/cosmology/tests/test_cosmology.py
func TestWCDMELcdm(t *testing.T) {
	cos := WCDM{Om0: 0.27, Ol0: 0.73, W0: -1, H0: 70, Tcmb0: 0.}

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

func TestWCDMDistanceModulus(t *testing.T) {
	cos := WCDM{Om0: 0.3, Ol0: 0.7, W0: -1.2, H0: 70, Tcmb0: 0.}
	//   wCDM(70, 0.3, 0.7, -1.2).distmod(z)
	exp_vec := answersWCDM["WCDMDistanceModulus"]
	runTests(cos.DistanceModulus, zWCDM, exp_vec, distTol, t)
}

func TestWCDMLuminosityDistanceFlatCDM(t *testing.T) {
	cos := WCDM{Om0: 0.3, Ol0: 0.7, W0: -1, H0: 70, Tcmb0: 0.}
	exp_vec := answersWCDM["WCDMLuminosityDistanceFlatCDM"]
	runTests(cos.LuminosityDistance, zWCDM, exp_vec, distTol, t)
}

func TestWCDMLuminosityDistanceFlat(t *testing.T) {
	cos := WCDM{Om0: 0.3, Ol0: 0.7, W0: -1.1, H0: 70, Tcmb0: 0.}
	//   wCDM(70, 0.3, 0.7, -1.1).luminosity_distance(z)
	exp_vec := answersWCDM["WCDMLuminosityDistanceFlat"]
	runTests(cos.LuminosityDistance, zWCDM, exp_vec, distTol, t)
}

func TestWCDMLuminosityDistanceNonflat(t *testing.T) {
	cos := WCDM{Om0: 0.3, Ol0: 0.6, W0: -0.8, H0: 70, Tcmb0: 0.}
	//   wCDM(70, 0.3, 0.6, -0.8).luminosity_distance(z)
	exp_vec := answersWCDM["WCDMLuminosityDistanceNonflat"]
	runTests(cos.LuminosityDistance, zWCDM, exp_vec, distTol, t)
}

func TestWCDMAngularDiameterDistance(t *testing.T) {
	cos := WCDM{Om0: 0.3, Ol0: 0.7, W0: -1, H0: 70, Tcmb0: 0.}
	exp_vec := answersWCDM["WCDMAngularDiameterDistance"]
	runTests(cos.AngularDiameterDistance, zWCDM, exp_vec, distTol, t)
}

func TestWCDMComovingDistanceNonflatOM(t *testing.T) {
	cos := WCDM{Om0: 0.3, Ol0: 0., W0: -1, H0: 70, Tcmb0: 0.}
	exp_vec := answersWCDM["WCDMComovingDistanceNonflatOM"]
	runTests(cos.ComovingDistance, zWCDM, exp_vec, distTol, t)
}

func TestWCDMComovingTransverseDistanceNonflatOM(t *testing.T) {
	cos := WCDM{Om0: 0.3, Ol0: 0., W0: -1, H0: 70, Tcmb0: 0.}
	exp_vec := answersWCDM["WCDMComovingTransverseDistanceNonflatOM"]
	runTests(cos.ComovingTransverseDistance, zWCDM, exp_vec, distTol, t)
}

func TestWCDMComovingDistanceEdS(t *testing.T) {
	cos := WCDM{Om0: 1.0, Ol0: 0., W0: -1, H0: 70, Tcmb0: 0.}
	exp_vec := answersWCDM["WCDMComovingDistanceEdS"]
	runTests(cos.ComovingDistance, zWCDM, exp_vec, distTol, t)
}

func TestWCDMComovingTransverseDistance(t *testing.T) {
	cos := WCDM{Om0: 0.3, Ol0: 0.7, W0: -1, H0: 70, Tcmb0: 0.}
	exp_vec := answersWCDM["WCDMComovingTransverseDistance"]
	runTests(cos.ComovingTransverseDistance, zWCDM, exp_vec, distTol, t)
}

func TestWCDMComovingDistanceZ1Z2Integrate(t *testing.T) {
	cos := WCDM{Om0: 0.3, Ol0: 0.7, W0: -1, H0: 70, Tcmb0: 0.}
	exp_vec := answersWCDM["WCDMComovingDistanceZ1Z2Integrate"]
	runTestsZ0Z2(cos.comovingDistanceZ1Z2Integrate, zWCDM, exp_vec, distTol, t)
}

func TestWCDMComovingDistanceZ1Z2Elliptic(t *testing.T) {
	cos := WCDM{Om0: 0.3, Ol0: 0.7, W0: -1, H0: 70, Tcmb0: 0.}
	exp_vec := answersWCDM["WCDMComovingDistanceZ1Z2Elliptic"]
	runTestsZ0Z2(cos.comovingDistanceZ1Z2Elliptic, zWCDM, exp_vec, distTol, t)
}

func TestWCDMLookbackTime(t *testing.T) {
	cos := WCDM{Om0: 0.3, Ol0: 0.7, W0: -1.2, H0: 70, Tcmb0: 0.}
	// Calculated via astropy.cosmology.wCDM(70, 0.3, 0.7, w0=-1.2).lookback_time
	exp_vec := answersWCDM["WCDMLookbackTime"]
	runTests(cos.LookbackTime, zWCDM, exp_vec, ageTol, t)
}

func TestWCDMLookbackTimeIntegrate(t *testing.T) {
	cos := WCDM{Om0: 0.3, Ol0: 0.7, W0: -1.1, H0: 70, Tcmb0: 0.}
	// Calculated via astropy.cosmology.wCDM(70, 0.3, 0.7, w0=-1.1).lookback_time
	exp_vec := answersWCDM["WCDMLookbackTimeIntegrate"]
	runTests(cos.lookbackTimeIntegrate, zWCDM, exp_vec, ageTol, t)
}

func TestWCDMLookbackTimeOM(t *testing.T) {
	cos := WCDM{Om0: 0.3, Ol0: 0., W0: -0.9, H0: 70, Tcmb0: 0.}
	// Calculated via astropy.cosmology.wCDM(70, 0.3, 0.7, w0=-0.9).lookback_time
	exp_vec := answersWCDM["WCDMLookbackTimeOM"]
	runTests(cos.LookbackTime, zWCDM, exp_vec, ageTol, t)
}

func TestWCDMAge(t *testing.T) {
	cos := WCDM{Om0: 0.3, Ol0: 0.6, W0: -1, H0: 70, Tcmb0: 0.}
	//   astropy.cosmology.WCDM(70, 0.3, 0.6).age(z)
	exp_vec := answersWCDM["WCDMAge"]
	runTests(cos.Age, zWCDM, exp_vec, ageTol, t)
}

func TestWCDMAgeFlatLCDM(t *testing.T) {
	cos := WCDM{Om0: 0.3, Ol0: 0.7, W0: -1, H0: 70, Tcmb0: 0.}
	//   astropy.cosmology.wCDM(70, 0.3).age(z)
	exp_vec := answersWCDM["WCDMAgeFlatLCDM"]
	runTests(cos.Age, zWCDM, exp_vec, ageTol, t)
}

func TestWCDMAgeIntegrate(t *testing.T) {
	cos := WCDM{Om0: 0.3, Ol0: 0.7, W0: -1, H0: 70, Tcmb0: 0.}
	//   astropy.cosmology.wCDM(70, 0.3).age(z)
	exp_vec := answersWCDM["WCDMAgeIntegrate"]
	runTests(cos.ageIntegrate, zWCDM, exp_vec, ageTol, t)
}

func TestWCDMAgeOM(t *testing.T) {
	cos := WCDM{Om0: 0.3, Ol0: 0., W0: -1, H0: 70, Tcmb0: 0.}
	//   astropy.cosmology.WCDM(70, 0.3, 0.).age(z)
	exp_vec := answersWCDM["WCDMAgeOM"]
	runTests(cos.Age, zWCDM, exp_vec, ageTol, t)
	runTests(cos.ageIntegrate, zWCDM, exp_vec, ageTol, t)
}

func TestWCDMAgeOL(t *testing.T) {
	cos := WCDM{Om0: 0, Ol0: 0.5, W0: -1, H0: 70, Tcmb0: 0.}
	//   astropy.cosmology.WCDM(70, 0.3, 0.).age(z)
	exp_vec := answersWCDM["WCDMAgeOL"]
	runTests(cos.Age, zWCDM, exp_vec, ageTol, t)
	runTests(cos.ageIntegrate, zWCDM, exp_vec, ageTol, t)
}

// Analytic case of Omega_Lambda = 0
func TestWCDMEOm(t *testing.T) {
	cos := WCDM{Om0: 1.0, Ol0: 0., W0: -1, H0: 70, Tcmb0: 0.}

	z_vec := []float64{1.0, 10.0, 500.0, 1000.0}
	exp_vec := make([]float64, len(z_vec))
	hubbleDistance := SpeedOfLightKmS / cos.H0
	for i, z := range z_vec {
		exp_vec[i] = 2.0 * hubbleDistance * (1 - math.Sqrt(1/(1+z)))
	}
	runTests(cos.ComovingDistance, z_vec, exp_vec, distTol, t)
}
