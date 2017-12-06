// Package cosmogo implements basic cosmology calculations in Go
//
// Provides E, Einv, LookbackTime, Age, LuminosityDistance, DistanceModulus,
// ComovingDistance, AngularDiameterDistance
//
// Equations and numerical formulae based on
//   Hogg, https://arxiv.org/abs/astro-ph/9905116
//   Feige, 1992, Astron. Nachr., 313, 139.
//   Eisenstein, 1997, https://arXiv.org/abs/astro-ph/9709054v2
//   Mészáros & Řípai 2013, A&A, 556, A13.
//   Baes, Camps, Van De Putte, 2017, MNRAS, 468, 927.
//   Kantowski, 2000, https://arxiv.org/abs/astro-ph/0002334
//   Thomas and Kantowski, 2000, PRD, 62, 103507.  Eq. 3
//
// Organizational thoughts based on code in astropy.cosmology
//   http://docs.astropy.org/en/stable/_modules/astropy/cosmology
package cosmo

const SpeedOfLightKmS = 299792.458 // km/s
// https://en.wikipedia.org/wiki/Parsec
//   Accessed 2017-12-01
const kmInAMpc = 3.08567758149137e19 // km/Mpc
//   365*24*3600 * one billion
const secInAGyr = 31557600 * 1e9 // s/Gyr

// Cosmology specifies all of the functions we need for cosmological calculation
type Cosmology interface {
	Age(z float64) (time float64)
	AgeFlatLCDM(z float64) (time float64)
	AgeIntegrate(z float64) (time float64)
	AgeOL(z float64) (time float64)
	AgeOM(z float64) (time float64)
	AngularDiameterDistance(z float64) (distance float64)
	ComovingDistance(z float64) (distance float64)
	ComovingDistanceOM(z float64) (distance float64)
	ComovingDistanceOMZ1Z2(z1, z2 float64) (distance float64)
	ComovingDistanceZ1Z2(z1, z2 float64) (distance float64)
	ComovingDistanceZ1Z2Elliptic(z1, z2 float64) (distance float64)
	ComovingDistanceZ1Z2Integrate(z1, z2 float64) (distance float64)
	ComovingTransverseDistance(z float64) (distance float64)
	ComovingTransverseDistanceZ1Z2(z1, z2 float64) (distance float64)
	DistanceModulus(z float64) (distmod float64)
	E(z float64) (ez float64)
	Einv(z float64) (invEz float64)
	HubbleDistance() float64
	LookbackTime(z float64) (time float64)
	LookbackTimeIntegrate(z float64) (time float64)
	LookbackTimeOL(z float64) (time float64)
	LookbackTimeOM(z float64) (time float64)
	LuminosityDistance(z float64) (distance float64)
}
