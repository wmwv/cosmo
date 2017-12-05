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

import (
	"gonum.org/v1/gonum/integrate/quad"
	"gonum.org/v1/gonum/mathext"
	"math"
)

const SpeedOfLightKmS = 299792.458 // km/s
// https://en.wikipedia.org/wiki/Parsec
//   Accessed 2017-12-01
const kmInAMpc = 3.08567758149137e19 // km/Mpc
//   365*24*3600 * one billion
const secInAGyr = 31557600 * 1e9 // s/Gyr

// Cosmology stores the key information needed for a given cosmology
type Cosmology struct {
	Om0     float64 // Matter Density at z=0
	Ol0     float64 // Vacuum Energy density Lambda at z=0
	Ok0     float64 // Curvature Density at z=0
	H0      float64 // Hubble constant at z=0.  [km/s/Mpc]
	W0      float64 // Dark energy equation-of-state parameter
	Ogamma0 float64 // Photon density
	Onu0    float64 // Neutrino density
	Tcmb0   float64 // Temperature of the CMB at z=0.  [K]
	//    nuToPhotonDensity float64 // Neutrino density / photon density
}

// DistanceModulus is the magnitude difference between 1 Mpc and
// the luminosity distance for the given z.
// z : redshift
//
// distmod : [mag]
func (cos *Cosmology) DistanceModulus(z float64) (distmod float64) {
	return 5*math.Log10(cos.LuminosityDistance(z)) + 25
}

// LuminosityDistance is the radius of effective sphere over which the light has spread out
// z : redshift
//
// distance : Mpc
func (cos *Cosmology) LuminosityDistance(z float64) (distance float64) {
	return (1 + z) * cos.ComovingTransverseDistance(z)
}

// AngularDistance is the ratio of physical transverse size to angular size
// z : redshift
//
// distance : Mpc/rad
func (cos *Cosmology) AngularDiameterDistance(z float64) (distance float64) {
	return cos.ComovingTransverseDistance(z) / (1 + z)
}

// ComovingTransverseDistance is the comoving distance at z as seen from z=0
// z : redshift
//
// distance : Mpc
func (cos *Cosmology) ComovingTransverseDistance(z float64) (distance float64) {
	return cos.ComovingTransverseDistanceZ1Z2(0, z)
}

// ComovingTransverseDistanceZ1Z2 is the comoving distance at z2 as seen from z1
// z : redshift
//
// distance : Mpc
func (cos *Cosmology) ComovingTransverseDistanceZ1Z2(z1, z2 float64) (distance float64) {
	comovingDistance := cos.ComovingDistanceZ1Z2(z1, z2)
	if cos.Ok0 == 0 {
		return comovingDistance
	}

	hubbleDistance := cos.HubbleDistance()
	return hubbleDistance /
		math.Sinh(math.Sqrt(cos.Ok0)*comovingDistance/hubbleDistance)
}

// HubbleDistance is the inverse of the Hubble parameter
//
// distance : Mpc
func (cos *Cosmology) HubbleDistance() float64 {
	return SpeedOfLightKmS / cos.H0
}

// ComovingDistance is the distance that is constant with the Hubble flow
// expressed in the physical distance at z=0.
// I.e., as the scale factor a = 1/(1+z) increases from 0.5 to 1,
// two objects separated by a proper distance of 10 Mpc at a=0.5 (z=1)
// will be separated by a proper distance of 10*2 Mpc at a=1.0 (z=0)
// The comoving distance between these objects is 20 Mpc
//
func (cos *Cosmology) ComovingDistance(z float64) (distance float64) {
	return cos.ComovingDistanceZ1Z2(0, z)
}

// ComovingDistanceZ1Z2Elliptic is the comoving distance between two z
// in a flat lambda CDM cosmology using elliptic integrals.
//
// See
//     Feige, 1992, Astron. Nachr., 313, 139.
//     Eisenstein, 1997, arXiv:9709054v2
//     Mészáros & Řípai 2013, A&A, 556, A13.
// and a useful summary in
//     Baes, Camps, Van De Putte, 2017, MNRAS, 468, 927.
func (cos *Cosmology) ComovingDistanceZ1Z2Elliptic(z1, z2 float64) (distance float64) {
	s := math.Pow((1-cos.Om0)/cos.Om0, 1./3)
	prefactor := (SpeedOfLightKmS / cos.H0) * (1 / math.Sqrt(s*cos.Om0))
	return prefactor * (TElliptic(s/(1+z1)) - TElliptic(s/(1+z2)))
}

// TElliptic uses elliptic integral of the first kind in Carlson form
//   to calculate the basic integral for cosmological distances
// gonum.org/v1/mathext/EllipticRF (Carlson form)
func TElliptic(s float64) float64 {
	m := (2 * math.Sqrt(s*s-s+1) / s) + (2 / s) - 1
	x := m
	y := m + 3 - 2*math.Sqrt(3)
	z := m + 3 + 2*math.Sqrt(3)
	return 4 * mathext.EllipticRF(x, y, z)
}

// ComovingDistanceZ1Z2Integrate is the comoving distance between two z
//   in a flat lambda CDM cosmology using fixed Gaussian quadrature integration.
func (cos *Cosmology) ComovingDistanceZ1Z2Integrate(z1, z2 float64) (distance float64) {
	n := 1000 // Integration will be n-point Gaussian quadrature
	return cos.HubbleDistance() * quad.Fixed(cos.Einv, z1, z2, n, nil, 0)
}

// ComovingDistanceZ1Z2 is the base function for calculation of comoving distances
//   Here is where the choice of fundamental calculation method is made:
//   Elliptic integral, quadrature integration, or analytic for special cases.
func (cos *Cosmology) ComovingDistanceZ1Z2(z1, z2 float64) (distance float64) {
	switch {
	// Test for Ol0==0 first so that (Om0, Ol0) = (1, 0)
	// is handled by the analytic solution
	// rather than the explicit integration.
	case cos.Ol0 == 0:
		return cos.ComovingDistanceOMZ1Z2(z1, z2)
	case cos.Om0 < 1:
		return cos.ComovingDistanceZ1Z2Elliptic(z1, z2)
	default:
		return cos.ComovingDistanceZ1Z2Integrate(z1, z2)
	}
}

// ComovingDistanceOM is the analytic case of Omega_total=Omega_M
func (cos *Cosmology) ComovingDistanceOM(z float64) (distance float64) {
	// Call the internal function that just takes direct arguments
	//   with nothing passed via the struct.
	return comovingDistanceOM(z, cos.Om0, cos.H0)
}

// ComovingDistanceOMZ1Z2 is the analytic case of Omega_total=Omega_M
//    for the distance between two redshifts.
//
// This *Z1Z2 form exists to parallel the other versions
//  and allow it to be a shortcut option in ComovingDistanceZ1Z2.
// Naively, it's twice as expensive to do this as (0, z2)
// But this is such a trivial calculation, it probably doesn't matter.
func (cos *Cosmology) ComovingDistanceOMZ1Z2(z1, z2 float64) (distance float64) {
	return comovingDistanceOM(z2, cos.Om0, cos.H0) -
		comovingDistanceOM(z1, cos.Om0, cos.H0)
}

// comovingDistanceOM is the analytic case of Omega_total=Omega_M
//
// Hogg, 1999
// Peebles, 1993
// Weinberg, 1972
// Mattig, 1958
// Transcribed from Kantowski 2000 (arXiv:0002334)
func comovingDistanceOM(z, Om0, H0 float64) (distance float64) {
	return (SpeedOfLightKmS / H0) *
		2 * (2 - Om0*(1-z) - (2-Om0)*math.Sqrt(1+Om0*z)) /
		((1 + z) * Om0 * Om0)
}

// LookbackTime is the time from redshift 0 to z.
//
// z : redshift
func (cos *Cosmology) LookbackTime(z float64) (time float64) {
	switch {
	case (cos.Ol0 == 0) && (0 < cos.Om0) && (cos.Om0 != 1):
		return cos.LookbackTimeOM(z)
	case (cos.Om0 == 0) && (0 < cos.Ol0) && (cos.Ol0 < 1):
		return cos.LookbackTimeOL(z)
	default:
		return cos.LookbackTimeIntegrate(z)
	}
}

// LookbackTimeIntegrate is the lookback time using explicit integration
//
// z : redshift
func (cos *Cosmology) LookbackTimeIntegrate(z float64) (time float64) {
	n := 1000 // Integration will be n-point Gaussian quadrature
	integrand := func(z float64) float64 { return cos.Einv(z) / (1 + z) }
	return hubbleTime(cos.H0) * quad.Fixed(integrand, 0, z, n, nil, 0)
}

// LookbackTimeOL is lookback time for dark-energy only Universe
//
// z : redshift
func (cos *Cosmology) LookbackTimeOL(z float64) (time float64) {
	return lookbackTimeOL(z, cos.Ol0, cos.H0)
}

// LookbackTimeOM is lookback time for matter only Universe
// All matter is non-relativistic.
//
// z : redshift
func (cos *Cosmology) LookbackTimeOM(z float64) (time float64) {
	return lookbackTimeOM(z, cos.Om0, cos.H0)
}

// lookbackTimeOL is lookback time for dark-energy + curvature only Universe
//
// z : redshift
// Ol0 : Omega_Lambda at z=0.
//       Dark energy density as a fraction of the critical density
// H0 : Hubble Parameter at z=0.  [km/s/Mpc]

// Equation is in many sources.  Sppecifically used
// Thomas and Kantowski, 2000, PRD, 62, 103507.  Eq. 3
func lookbackTimeOL(z, Ol0, H0 float64) (time float64) {
	return ageOL(0, Ol0, H0) - ageOL(z, Ol0, H0)
}

// lookbackTimeOM is lookback time for matter only + curvature Universe
// All matter is non-relativistic.
//
// z : redshift
// Om0 : Omega_M at z=0.
//       Matter density as a fraction of the critical density
//       All matter non-relatisvistic.
// H0 : Hubble Parameter at z=0.  [km/s/Mpc]
//
// Equation is in many sources.  Sppecifically used
// Thomas and Kantowski, 2000, PRD, 62, 103507.  Eq. 2
func lookbackTimeOM(z, Om0, H0 float64) (time float64) {
	return ageOM(0, Om0, H0) - ageOM(z, Om0, H0)
}

// Age is the time from redshift ∞ to z.
//
// z : redshift
//
// Method of Cosmology.  Requires defined Ol0, Om0, H0.
func (cos *Cosmology) Age(z float64) (time float64) {
	switch {
	case cos.Om0+cos.Ol0 == 1:
		return cos.AgeFlatLCDM(z)
	case (cos.Ol0 == 0) && (0 < cos.Om0) && (cos.Om0 != 1):
		return cos.AgeOM(z)
	case (cos.Om0 == 0) && (0 < cos.Ol0) && (cos.Ol0 < 1):
		return cos.AgeOL(z)
	default:
		return cos.AgeIntegrate(z)
	}
}

// AgeFlatLCDM is the time from redshift ∞ to z
// in a flat LCDM cosmology.
//
// Equation is in many sources.
// I took this from Thomas and Kantowski, 2000 PRD, 62, 103507.
func (cos *Cosmology) AgeFlatLCDM(z float64) (time float64) {
	return hubbleTime(cos.H0) * 2. / 3 / math.Sqrt(1-cos.Om0) *
		math.Asinh(math.Sqrt((1/cos.Om0-1)/math.Pow(1+z, 3)))
}

// AgeIntegrate is the time from redshift ∞ to z
// using explicit integration.
//
// Basic integrand can be found in many texts
// I happened to copy this from
// Thomas and Kantowski, 2000, PRD, 62, 103507.  Eq. 1.
// Current implementation is fixed quadrature using mathext.integrate.quad.Fixed
func (cos *Cosmology) AgeIntegrate(z float64) (time float64) {
	n := 1000 // Integration will be n-point Gaussian quadrature
	integrand := func(z float64) float64 {
		denom := (1 + z) * math.Sqrt((1+z)*(1+z)*(1+cos.Om0*z)-z*(2+z)*cos.Ol0)
		return 1 / denom
	}
	// When given math.Inf(), quad.Fixed automatically redefines variables
	// to successfully do the numerical integration.
	return hubbleTime(cos.H0) * quad.Fixed(integrand, z, math.Inf(1), n, nil, 0)
}

// AgeOL is the time from redshift ∞ to z
// with only constant dark energy and curvature.
// z : redshift
func (cos *Cosmology) AgeOL(z float64) (time float64) {
	return ageOL(z, cos.Ol0, cos.H0)
}

// AgeOL is the time from redshift ∞ to z
// with only non-relativistic matter and curvature.
// z : redshift
func (cos *Cosmology) AgeOM(z float64) (time float64) {
	return ageOM(z, cos.Om0, cos.H0)
}

// Calculate the Hubble time, c/H0.
//
// H0 : Hubble parameter at z=0.  [km/s/Mpc]
// Returns time in Gyr
func hubbleTime(H0 float64) (time float64) {
	hubbleTime := (1 / H0)  // 1/(km/s/Mpc) = Mpc s / km
	hubbleTime *= kmInAMpc  // s
	hubbleTime /= secInAGyr // Gyr

	return hubbleTime
}

// ageOL is the time from redshift ∞ to z
// with only constant dark energy and curvature.
// Bare function version.  Not method of struct Cosmology, just takes 3 floats.
// z : redshift
// Ol0 : Omega_Lambda at z=0.
//       Dark energy density as a fraction of the critical density
// H0 : Hubble Parameter at z=0.  [km/s/Mpc]
//
// Equation is in many sources.  Sppecifically used
// Thomas and Kantowski, 2000, PRD, 62, 103507.  Eq. 3
func ageOL(z, Ol0, H0 float64) (time float64) {
	return hubbleTime(H0) * (1 / math.Sqrt(Ol0)) *
		math.Asinh(1/((1+z)*math.Sqrt((1/Ol0)-1)))
}

// AgeOM is the time from redshift ∞ to z
// with only non-relativisitc matter and curvature.
// Bare function version.  Not method of struct Cosmology, just takes 3 floats.
// z : redshift
// Om0 : Omega_M at z=0.
//       Matter density as a fraction of the critical density
//       All matter non-relatisvistic.
// H0 : Hubble Parameter at z=0.  [km/s/Mpc]
//
// Equation is in many sources.  Specifically used
// Thomas and Kantowski, 2000, PRD, 62, 103507.  Eq. 2
func ageOM(z, Om0, H0 float64) (time float64) {
	return hubbleTime(H0) *
		(math.Sqrt(1+Om0*z)/((1-Om0)*(1+z)) -
			Om0*math.Pow(1-Om0, -3./2)*math.Asinh(math.Sqrt((1/Om0-1)/(1+z))))
}

// E is the Hubble parameter as a fraction of its present value.
// E.g., Hogg arXiv:9905116  Eq. 14
func (cos *Cosmology) E(z float64) (ez float64) {
	oR := cos.Ogamma0 + cos.Onu0
	var deScale float64
	// TODO
	// Consider an if or switch on the value of cos.W0
	// Do performance testing to see in what circumstances it matters.
	switch cos.W0 {
	case -1:
		deScale = 1
	default:
		deScale = math.Pow(1+z, 3*(1+cos.W0))
	}
	ez = math.Sqrt((1+z)*(1+z)*((oR*(1+z)+cos.Om0)*(1+z)+cos.Ok0) + cos.Ol0*deScale)
	return ez
}

// Einv is the inverse Hubble parameter
// Implementation is just to return E(z)
func (cos *Cosmology) Einv(z float64) (invEz float64) {
	// 1/Sqrt() is not notably slower than Pow(-0.5)
	//
	// Pow(-0.5) is in fact implemented as 1/Sqrt() in math.pow.go
	// func pow(x, y float64) float64 {
	//    [...]
	// case y == -0.5:
	//    return 1 / Sqrt(x)
	//
	// Thus we just return the inverse of E(z) instead of rewriting out here.
	return 1 / cos.E(z)
}

// Evec is vectorized form of 'E'.
// I haven't figured out whether this is useful or makes sense
// in a Go framework, which looping over functions is more expected.
// thank in IDL, Matlab, or Python numpy+scipy worlds.
func (cos *Cosmology) Evec(z []float64) (ez []float64) {
	for _, z := range z {
		ez = append(ez, cos.E(z))
	}
	return ez
}
