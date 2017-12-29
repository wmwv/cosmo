package cosmo

import (
	"gonum.org/v1/gonum/mathext"
	"math"
)

// lookbackTimeOL is lookback time for dark-energy + curvature only Universe
//   z : redshift
//   Ol0 : Omega_Lambda at z=0.
//         Dark energy density as a fraction of the critical density
//   H0 : Hubble Parameter at z=0.  [km/s/Mpc]
func lookbackTimeOL(z, Ol0, H0 float64) (timeGyr float64) {
	return ageOL(0, Ol0, H0) - ageOL(z, Ol0, H0)
}

// lookbackTimeOM is lookback time for matter only + curvature Universe
// All matter is non-relativistic.
//   z : redshift
//   Om0 : Omega_M at z=0.
//         Matter density as a fraction of the critical density
//         All matter non-relatisvistic.
//   H0 : Hubble Parameter at z=0.  [km/s/Mpc]
func lookbackTimeOM(z, Om0, H0 float64) (timeGyr float64) {
	return ageOM(0, Om0, H0) - ageOM(z, Om0, H0)
}

// Calculate the Hubble time, c/H0.
//   H0 : Hubble parameter at z=0.  [km/s/Mpc]
func hubbleTime(H0 float64) (timeGyr float64) {
	hubbleTime := (1 / H0)  // 1/(km/s/Mpc) = Mpc s / km
	hubbleTime *= kmInAMpc  // s
	hubbleTime /= secInAGyr // Gyr

	return hubbleTime
}

// ageOL is the time from redshift ∞ to z
// with only constant dark energy and curvature.
// Bare function version.  Not method of struct LambdaCDM, just takes 3 floats.
//   z : redshift
//   Ol0 : Omega_Lambda at z=0.
//         Dark energy density as a fraction of the critical density
//   H0 : Hubble Parameter at z=0.  [km/s/Mpc]
//
// Equation is in many sources.  Sppecifically used
// Thomas and Kantowski, 2000, PRD, 62, 103507.  Eq. 3
func ageOL(z, Ol0, H0 float64) (timeGyr float64) {
	return hubbleTime(H0) * (1 / math.Sqrt(Ol0)) *
		math.Asinh(1/((1+z)*math.Sqrt((1/Ol0)-1)))
}

// AgeOM is the time from redshift ∞ to z
// with only non-relativisitc matter and curvature.
// Bare function version.  Not method of struct LambdaCDM, just takes 3 floats.
//   z : redshift
//   Om0 : Omega_M at z=0.
//         Matter density as a fraction of the critical density
//         All matter non-relatisvistic.
//   H0 : Hubble Parameter at z=0.  [km/s/Mpc]
//
// Equation is in many sources.  Specifically used
// Thomas and Kantowski, 2000, PRD, 62, 103507.  Eq. 2
func ageOM(z, Om0, H0 float64) (timeGyr float64) {
	if Om0 == 1 {
		return (2. / 3) * hubbleTime(H0) * math.Pow(1+z, -3./2)
	}
	return hubbleTime(H0) *
		(math.Sqrt(1+Om0*z)/((1-Om0)*(1+z)) -
			Om0*math.Pow(1-Om0, -3./2)*math.Asinh(math.Sqrt((1/Om0-1)/(1+z))))
}

// ageFlatLCDM is the time from redshift ∞ to z
// with only non-relativistic matter and dark energy.  No curvature: Om0+Ol0=1
//
//   z : redshift
//   Om0 : Omega_M at z=0.
//         Matter density as a fraction of the critical density
//         All matter non-relatisvistic.
//   H0 : Hubble Parameter at z=0.  [km/s/Mpc]
//
// Equation is in many sources.  Specifically used
// Thomas and Kantowski, 2000, PRD, 62, 103507.
func ageFlatLCDM(z, Om0, H0 float64) (timeGyr float64) {
	if Om0 == 1 {
		return (2. / 3) * hubbleTime(H0) * math.Pow(1+z, -3./2)
	}
	return hubbleTime(H0) * 2. / 3 / math.Sqrt(1-Om0) *
		math.Asinh(math.Sqrt((1/Om0-1)/math.Pow(1+z, 3)))
}

// comovingTransverseDistanceOM is the case of Omega_M+Omega_K=1
//
//   z : redshift
//   Om0 : Omega_M at z=0.
//         Matter density as a fraction of the critical density
//         All matter non-relatisvistic.
//   H0 : Hubble Parameter at z=0.  [km/s/Mpc]
//
// Hogg, arXiv:9905116
// Peebles, 1993
// Weinberg, 1972
// Mattig, 1958
// Transcribed from Kantowski 2000 (arXiv:0002334)
func comovingTransverseDistanceOM(z, Om0, H0 float64) (distanceMpc float64) {
	return (SpeedOfLightKmS / H0) *
		2 * (2 - Om0*(1-z) - (2-Om0)*math.Sqrt(1+Om0*z)) /
		((1 + z) * Om0 * Om0)
}

// comovingDistanceOM is the case of Omega_M+Omega_K=1
//
// If Omega_K=0, then comovingDistance == comovingTransverseDistance
func comovingDistanceOM(z, Om0, H0 float64) (distanceMpc float64) {
	comovingTransverseDistance := comovingTransverseDistanceOM(z, Om0, H0)
	Ok0 := 1 - Om0
	if Ok0 == 0 {
		return comovingTransverseDistance
	}

	hubbleDistance := SpeedOfLightKmS / H0
	hdk := hubbleDistance / math.Sqrt(Ok0)
	return hdk * math.Asinh(comovingTransverseDistance/hdk)
}

// tElliptic uses elliptic integral of the first kind in Carlson form
//   to calculate the basic integral for cosmological distances
// gonum.org/v1/mathext/EllipticRF (Carlson form)
func tElliptic(s float64) float64 {
	m := (2 * math.Sqrt(s*s-s+1) / s) + (2 / s) - 1
	x := m
	y := m + 3 - 2*math.Sqrt(3)
	z := m + 3 + 2*math.Sqrt(3)
	return 4 * mathext.EllipticRF(x, y, z)
}
