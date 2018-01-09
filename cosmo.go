// Package cosmo implements basic cosmology calculations in Go.
//
// FLRW is the basic interface type that defines the key cosmological functions
// common to all supported cosmologies.
//
// The Friedmann-Lemaître-Robertson-Walker (FLRW) metric is the general form that
// all homogenous, isotropic, connected cosmologies follow.
//
// https://en.wikipedia.org/wiki/Friedmann-Lema%C3%AEtre-Robertson-Walker_metric
//
// Provides:
//   FlatLCDM  (H0, OM); OL = 1-OM, OK=0
//   LambdaCDM (H0, OM, OL, OK); w = -1
//   WCDM      (H0, OM, OL, W); w = w0
//   WACDM     (H0, OM, OL, W0, WA); w = w0 + w_a * (1-a)
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
// Performance Notes:
//
// 1. Calculating a single line-of-sight comoving distance takes from 890ns - 260µs,
// depending on the complexity of the cosmology.
// Analytic cases take ~1µs, while explicit integration is ~200µs.
//   FlatLCDM      892ns  (analytic for OM<1)
//   LambdaCDM     139µs
//   WCDM          250µs
//   WACDM         261µs
// These numbers are based on the output of `go test -bench ComovingDistance`
// run on a 2015 MacBook Air: dual-core 2.2 GHz Intel Core i7, 8 GB 1600 MHz DDR3;
// with Mac OS X 10.13.2 and go 1.9.2
//
// 2. The current types FlatLCDM, LambdaCDM, WCDM, WACDM
// implement their methods as value receivers.
// There's a mild performance hit for using value receivers instead of pointer receivers.
// This performance penalty is 40% for individual calls to E or Einv
// but this penalty is only 1-2% for calls to the general *Distance methods.
//
// For now the use case model seems more amenable to value receivers
// -- conceptually, the cosmologies should be immutable --
// and the performance penalty is acceptable.
package cosmo

// SpeedOfLightKmS is speed of light in kilometers/second
// Useful in calculating Hubble distance (c/H0),
// which is the basic prefactor for any distance measure.
const SpeedOfLightKmS = 299792.458 // km/s
// https://en.wikipedia.org/wiki/Parsec
//   Accessed 2017-12-01
const kmInAMpc = 3.08567758149137e19 // km/Mpc
//   365*24*3600 * one billion
const secInAGyr = 31557600 * 1e9 // s/Gyr

// FLRW specifies the cosmological calculations to be available from
// Friedmann-Lemaître-Robertson-Walker metrics.
type FLRW interface {
	Age(z float64) (timeGyr float64)
	AngularDiameterDistance(z float64) (distanceMpc float64)
	ComovingDistance(z float64) (distanceMpc float64)
	ComovingDistanceZ1Z2(z1, z2 float64) (distanceMpc float64)
	ComovingTransverseDistance(z float64) (distanceMpc float64)
	ComovingTransverseDistanceZ1Z2(z1, z2 float64) (distanceMpc float64)
	DistanceModulus(z float64) (distanceModulusMag float64)
	E(z float64) (fractionalHubbleParameter float64)
	Einv(z float64) (invFractionalHubbleParameter float64)
	HubbleDistance() (distanceMpc float64)
	LookbackTime(z float64) (timeGyr float64)
	LuminosityDistance(z float64) (distanceMpc float64)
	Ok0() (curvatureDensity float64)
}
