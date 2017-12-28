// Package cosmo implements basic cosmology calculations in Go
//
// FLRW is the basic interface type that defines the functions to be supported
// The Friedmann-Lemaître-Robertson-Walker metric is the general form that
// all homogenous, isotropic, connected cosmologies follow.
// https://en.wikipedia.org/wiki/Friedmann-Lema%C3%AEtre-Robertson-Walker_metric
//
// Provides E, Einv, LookbackTime, Age, LuminosityDistance, DistanceModulus,
// ComovingDistance, ComovingTransverseDistance, AngularDiameterDistance
// LookbackTime, HubbleDistance
//
//   FlatLCDM  (OM, OL, OK) = (OM, 1-OM, 0)
//   LambdaCDM (OM, OL, OK); w = -1
//   WCDM      (OM, OL, W); w = w0
//   WACDM     (OM, OL, W0, WA); w = w0 + w_a * (1-a)
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
//
// Performance Note:
//
// 1. The current types FlatLCDM, LambdaCDM, WCDM, WACDM
// implement their methods as value receivers.
// There's a mild performance hit for using value receivers instead of pointer receivers.
// This performance penalty is 40% for individual calls to E or Einv
// but this penalty is only 1-2% for calls to the general *Distance methods.
//
// For now the use case model seems more amenable to value receivers
// and the performance penalty is acceptable.
package cosmo

const SpeedOfLightKmS = 299792.458 // km/s
// https://en.wikipedia.org/wiki/Parsec
//   Accessed 2017-12-01
const kmInAMpc = 3.08567758149137e19 // km/Mpc
//   365*24*3600 * one billion
const secInAGyr = 31557600 * 1e9 // s/Gyr

// FLRW specifies the functions cosmological calculations in
// Friedmann-Lemaître-Robertson-Walker metrics.
type FLRW interface {
	Age(z float64) (time float64)
	AngularDiameterDistance(z float64) (distance float64)
	ComovingDistance(z float64) (distance float64)
	ComovingDistanceZ1Z2(z1, z2 float64) (distance float64)
	ComovingTransverseDistance(z float64) (distance float64)
	ComovingTransverseDistanceZ1Z2(z1, z2 float64) (distance float64)
	DistanceModulus(z float64) (distmod float64)
	E(z float64) (ez float64)
	Einv(z float64) (invEz float64)
	HubbleDistance() (distance float64)
	LookbackTime(z float64) (time float64)
	LookbackTimeIntegrate(z float64) (time float64)
	LuminosityDistance(z float64) (distance float64)
}
