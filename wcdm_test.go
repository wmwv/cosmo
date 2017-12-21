package cosmo

import (
	"fmt"
	"gonum.org/v1/gonum/floats"
	"math"
	"runtime"
	"testing"
)

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
	exp := 1.7
	obs := cos.E(1.0)
	tol := 1e-9
	if !floats.EqualWithinAbs(obs, exp, tol) {
		t.Errorf("Failed flat wCDM E(z) test.  Expected %f, return %f",
			exp, obs)
	}

	exp = 1 / 1.7
	obs = cos.Einv(1.0)
	if !floats.EqualWithinAbs(obs, exp, tol) {
		t.Errorf("Failed flat wCDM Einv(z) test.  Expected %f, return %f",
			exp, obs)
	}
}

func TestWCDMDistanceModulus(t *testing.T) {
	var obs float64
	cos := WCDM{Om0: 0.3, Ol0: 0.7, W0: -1.2, H0: 70, Tcmb0: 0.}

	tol := 1e-8
	//  z_vec := []float64{0.2, 0.4, 0.9, 1.2}
	//  exp_vec := []float64{971.667, 2141.67, 5685.96, 8107.41}
	// Calculated via
	//   from astropy.cosmology import wCDM
	//   z = np.asarray([0.5, 1.0, 2.0, 3.0])
	//   wCDM(70, 0.3, 0.7, -1.2).distmod(z)

	z_vec := []float64{0.5, 1.0, 2.0, 3.0}
	exp_vec := []float64{42.32710911, 44.17957201, 46.03118144, 47.09228735}
	for i, z := range z_vec {
		obs = cos.DistanceModulus(z)
		if !floats.EqualWithinAbs(obs, exp_vec[i], tol) {
			t.Errorf("Failed flat wCDM luminosity distance test."+
				"  Expected %f, return %f",
				exp_vec[i], obs)
		}

	}
}

func TestWCDMLuminosityDistanceFlatCDM(t *testing.T) {
	var obs float64
	cos := WCDM{Om0: 0.3, Ol0: 0.7, W0: -1, H0: 70, Tcmb0: 0.}

	tol := 1e-6
	z_vec := []float64{0.5, 1.0, 2.0, 3.0}
	exp_vec := []float64{2832.9380939, 6607.65761177, 15539.58622323, 25422.74174519}
	for i, z := range z_vec {
		obs = cos.LuminosityDistance(z)
		if !floats.EqualWithinAbs(obs, exp_vec[i], tol) {
			t.Errorf("Failed flat wCDM luminosity distance test."+
				"  Expected %f, return %f",
				exp_vec[i], obs)
		}

	}
}

func TestWCDMLuminosityDistanceFlat(t *testing.T) {
	var obs float64
	cos := WCDM{Om0: 0.3, Ol0: 0.7, W0: -1.1, H0: 70, Tcmb0: 0.}

	tol := 1e-6
	// Calculated via
	//   from astropy.cosmology import LambdaCDM
	//   z = np.asarray([0.5, 1.0, 2.0, 3.0])
	//   wCDM(70, 0.3, 0.7, -1.1).luminosity_distance(z)
	z_vec := []float64{0.5, 1.0, 2.0, 3.0}
	exp_vec := []float64{2877.10314183, 6734.38177991, 15823.59621899, 25841.56448508}
	for i, z := range z_vec {
		obs = cos.LuminosityDistance(z)
		if !floats.EqualWithinAbs(obs, exp_vec[i], tol) {
			t.Errorf("Failed nonLambda, flat wCDM luminosity distance test."+
				"  Expected %f, return %f",
				exp_vec[i], obs)
		}

	}
}

func test(cos_func func(float64) float64, z, exp, tol float64, t *testing.T) {
	var test_description, test_line string

	pc, file, no, ok := runtime.Caller(1)
	if ok {
		details := runtime.FuncForPC(pc)
		test_description = details.Name()
		test_line = fmt.Sprintf("%s#%d", file, no)
	} else {
		test_description = "CAN'T DETERMINE TEST NAME"
		test_line = "CAN'T DETERMINE TEST LINE"
	}
	obs := cos_func(z)
	if !floats.EqualWithinAbs(obs, exp, tol) {
		t.Errorf("Failed %s at\n %s\n"+"  Expected %f, return %f",
			test_description, test_line, exp, obs)
	}
}

func TestWCDMLuminosityDistanceNonflat(t *testing.T) {
	cos := WCDM{Om0: 0.3, Ol0: 0.6, W0: -0.8, H0: 70, Tcmb0: 0.}

	tol := 1e-6
	// Calculated via
	//   from astropy.cosmology import LambdaCDM
	//   z = np.asarray([0.5, 1.0, 2.0, 3.0])
	//   wCDM(70, 0.3, 0.6, -0.8).luminosity_distance(z)
	z_vec := []float64{0.5, 1.0, 2.0, 3.0}
	exp_vec := []float64{2713.4660301, 6257.24866642, 14794.59911147, 24496.30592953}
	for i, z := range z_vec {
		test(cos.LuminosityDistance, z, exp_vec[i], tol, t)
	}
}

func TestWCDMAngularDiameterDistance(t *testing.T) {
	cos := WCDM{Om0: 0.3, Ol0: 0.7, W0: -1, H0: 70, Tcmb0: 0.}

	tol := 1e-6
	z_vec := []float64{0.5, 1.0, 2.0, 3.0}
	exp_vec := []float64{1259.08359729, 1651.91440294, 1726.62069147, 1588.92135907}
	for i, z := range z_vec {
		test(cos.AngularDiameterDistance, z, exp_vec[i], tol, t)
	}
}

func TestWCDMComovingTransverseDistance(t *testing.T) {
	cos := WCDM{Om0: 0.3, Ol0: 0.7, W0: -1, H0: 70, Tcmb0: 0.}

	tol := 1e-6
	//  z_vec := []float64{0.2, 0.4, 0.9, 1.2}
	//  exp_vec := []float64{971.667, 2141.67, 5685.96, 8107.41}
	z_vec := []float64{0.5, 1.0, 2.0, 3.0}
	exp_vec := []float64{1888.62539593, 3303.82880589, 5179.86207441, 6355.6854363}
	for i, z := range z_vec {
		test(cos.ComovingTransverseDistance, z, exp_vec[i], tol, t)
	}
}

func TestWCDMComovingDistanceZ1Z2Integrate(t *testing.T) {
	var obs float64
	cos := WCDM{Om0: 0.3, Ol0: 0.7, W0: -1, H0: 70, Tcmb0: 0.}

	tol := 1e-6
	//  z_vec := []float64{0.2, 0.4, 0.9, 1.2}
	//  exp_vec := []float64{971.667, 2141.67, 5685.96, 8107.41}
	z_vec := []float64{0.5, 1.0, 2.0, 3.0}
	exp_vec := []float64{1888.62539593, 3303.82880589, 5179.86207441, 6355.6854363}
	for i, z := range z_vec {
		obs = cos.ComovingDistanceZ1Z2Integrate(0, z)
		if !floats.EqualWithinAbs(obs, exp_vec[i], tol) {
			t.Errorf("Failed flat wCDM comoving distance elliptic test."+
				"  Expected %f, return %f",
				exp_vec[i], obs)
		}

	}
}

func TestWCDMComovingDistanceZ1Z2Elliptic(t *testing.T) {
	var obs float64
	cos := WCDM{Om0: 0.3, Ol0: 0.7, W0: -1, H0: 70, Tcmb0: 0.}

	tol := 1e-6
	//  z_vec := []float64{0.2, 0.4, 0.9, 1.2}
	//  exp_vec := []float64{971.667, 2141.67, 5685.96, 8107.41}
	z_vec := []float64{0.5, 1.0, 2.0, 3.0}
	exp_vec := []float64{1888.62539593, 3303.82880589, 5179.86207441, 6355.6854363}
	for i, z := range z_vec {
		obs = cos.ComovingDistanceZ1Z2Elliptic(0, z)
		if !floats.EqualWithinAbs(obs, exp_vec[i], tol) {
			t.Errorf("Failed flat wCDM comoving distance elliptic test."+
				"  Expected %f, return %f",
				exp_vec[i], obs)
		}

	}
}

func TestWCDMLookbackTime(t *testing.T) {
	cos := WCDM{Om0: 0.3, Ol0: 0.7, W0: -1.2, H0: 70, Tcmb0: 0.}

	tol := 1e-6
	z_vec := []float64{0.5, 1.0, 2.0, 3.0}
	// Calculated via astropy.cosmology.wCDM(70, 0.3, 0.7, w0=-1.2).lookback_time
	exp_vec := []float64{5.18796426, 7.98542226, 10.58842012, 11.71902479}
	for i, z := range z_vec {
		test(cos.LookbackTime, z, exp_vec[i], tol, t)
	}
}

func TestWCDMLookbackTimeIntegrate(t *testing.T) {
	cos := WCDM{Om0: 0.3, Ol0: 0.7, W0: -1.1, H0: 70, Tcmb0: 0.}

	tol := 1e-6
	z_vec := []float64{0.5, 1.0, 2.0, 3.0}
	// Calculated via astropy.cosmology.wCDM(70, 0.3, 0.7, w0=-1.1).lookback_time
	exp_vec := []float64{5.11509518, 7.85406053, 10.42213038, 11.54588106}
	for i, z := range z_vec {
		test(cos.LookbackTimeIntegrate, z, exp_vec[i], tol, t)
	}
}

func TestWCDMLookbackTimeOM(t *testing.T) {
	cos := WCDM{Om0: 0.3, Ol0: 0., W0: -0.9, H0: 70, Tcmb0: 0.}

	tol := 1e-6
	z_vec := []float64{0.5, 1.0, 2.0, 3.0}
	// Calculated via astropy.cosmology.wCDM(70, 0.3, 0.7, w0=-0.9).lookback_time
	exp_vec := []float64{4.51471693, 6.62532254, 8.57486509, 9.45923582}
	for i, z := range z_vec {
		test(cos.LookbackTimeOM, z, exp_vec[i], tol, t)
	}
}

func TestWCDMLookbackTimeOL(t *testing.T) {
	cos := WCDM{Om0: 0., Ol0: 0.5, W0: -1, H0: 70, Tcmb0: 0.}

	tol := 1e-6
	z_vec := []float64{0.5, 1.0, 2.0, 3.0}
	// Calculated via astropy.cosmology.wCDM(70, 0.3).lookback_time
	exp_vec := []float64{5.0616361, 7.90494991, 10.94241739, 12.52244605}
	for i, z := range z_vec {
		test(cos.LookbackTimeOL, z, exp_vec[i], tol, t)
	}
}

func TestWCDMAge(t *testing.T) {
	cos := WCDM{Om0: 0.3, Ol0: 0.6, W0: -1, H0: 70, Tcmb0: 0.}

	tol := 1e-6

	z_vec := []float64{0.5, 1.0, 2.0, 3.0}
	// Calculated via
	//   import astropy.cosmology
	//   z = [0.5, 1.0, 2.0, 3.0]
	//   astropy.cosmology.WCDM(70, 0.3, 0.6).age(z)
	exp_vec := []float64{8.11137578, 5.54558439, 3.13456008, 2.06445301}

	for i, z := range z_vec {
		test(cos.Age, z, exp_vec[i], tol, t)
	}
}

func TestWCDMAgeFlatLCDM(t *testing.T) {
	cos := WCDM{Om0: 0.3, Ol0: 0.7, W0: -1, H0: 70, Tcmb0: 0.}

	tol := 1e-6
	z_vec := []float64{0.5, 1.0, 2.0, 3.0}

	// Calculated via
	//   import astropy.cosmology
	//   z = [0.5, 1.0, 2.0, 3.0]
	//   astropy.cosmology.wCDM(70, 0.3).age(z)
	exp_vec := []float64{8.42634602, 5.75164694, 3.22662706, 2.11252719}

	for i, z := range z_vec {
		test(cos.Age, z, exp_vec[i], tol, t)
	}
}

func TestWCDMAgeIntegrate(t *testing.T) {
	cos := WCDM{Om0: 0.3, Ol0: 0.7, W0: -1, H0: 70, Tcmb0: 0.}

	tol := 1e-6
	z_vec := []float64{0.5, 1.0, 2.0, 3.0}

	// Calculated via
	//   import astropy.cosmology
	//   z = [0.5, 1.0, 2.0, 3.0]
	//   astropy.cosmology.wCDM(70, 0.3).age(z)
	exp_vec := []float64{8.42634602, 5.75164694, 3.22662706, 2.11252719}

	for i, z := range z_vec {
		test(cos.AgeIntegrate, z, exp_vec[i], tol, t)
	}
}

func TestWCDMAgeOM(t *testing.T) {
	cos := WCDM{Om0: 0.3, Ol0: 0., W0: -1, H0: 70, Tcmb0: 0.}

	tol := 1e-6
	z_vec := []float64{0.5, 1.0, 2.0, 3.0}
	// Calculated via
	//   import astropy.cosmology
	//   z = [0.5, 1.0, 2.0, 3.0]
	//   astropy.cosmology.WCDM(70, 0.3, 0.).age(z)
	exp_vec := []float64{6.78287955, 4.67227393, 2.72273139, 1.83836065}
	for i, z := range z_vec {
		test(cos.AgeOM, z, exp_vec[i], tol, t)
		test(cos.AgeIntegrate, z, exp_vec[i], tol, t)
	}
}

// Analytic case of Omega_Lambda = 0
func TestWCDMEOm(t *testing.T) {
	var exp float64
	cos := WCDM{Om0: 1.0, Ol0: 0., W0: -1, H0: 70, Tcmb0: 0.}

	tol := 1e-9
	z_vec := []float64{1.0, 10.0, 500.0, 1000.0}
	hubbleDistance := SpeedOfLightKmS / cos.H0
	for _, z := range z_vec {
		exp = 2.0 * hubbleDistance * (1 - math.Sqrt(1/(1+z)))
		test(cos.ComovingDistance, z, exp, tol, t)
	}
}

func TestWCDMEvecLcdm(t *testing.T) {
	cos := WCDM{Om0: 0.27, Ol0: 0.73, W0: -1, H0: 70, Tcmb0: 0.}
	// Check array
	z := []float64{0.5, 1.0}
	// Calculated using astropy.cosmology.wCDM (v1.3.2)
	exp := []float64{1.2811127975318957, 1.7}
	tol := 1e-9
	obs := cos.Evec(z)
	if !floats.EqualApprox(obs, exp, tol) {
		t.Errorf("Failed array float LCDM test.  Expected %v, return %v", exp, obs)
	}

}
