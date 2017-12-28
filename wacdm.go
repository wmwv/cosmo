package cosmo

import (
	"gonum.org/v1/gonum/integrate/quad"
	"math"
)

// WACDM provides cosmological distances, age, and look-back time
// for a w(a) cosmology:
// matter, dark energy, and curvature,
// with a
// w = w0 + wa * (1-a)
// equation-of-state parameter for dark energy.
//
// No radiation pressure or neutrino contributions.
type WACDM struct {
	Om0     float64 // Matter Density at z=0
	Ol0     float64 // Dark Energy density Lambda at z=0
	W0      float64 // Dark energy equation-of-state parameter, w0 + wa*(1-a) = p/rho
	WA      float64 // Dark energy equation-of-state parameter, w0 + wa*(1-a) = p/rho
	H0      float64 // Hubble constant at z=0.  [km/s/Mpc]
	Ogamma0 float64 // Photon density
	Onu0    float64 // Neutrino density
	Tcmb0   float64 // Temperature of the CMB at z=0.  [K]
	//    nuToPhotonDensity float64 // Neutrino density / photon density
}

// DistanceModulus is the magnitude difference between 1 Mpc and
// the luminosity distance for the given z.
//   z : redshift
//   distmod : distance modulus [mag]
func (cos WACDM) DistanceModulus(z float64) (distmod float64) {
	return 5*math.Log10(cos.LuminosityDistance(z)) + 25
}

// LuminosityDistance is the radius of effective sphere over which the light has spread out
//   z : redshift
//   distance : [Mpc]
func (cos WACDM) LuminosityDistance(z float64) (distance float64) {
	return (1 + z) * cos.ComovingTransverseDistance(z)
}

// AngularDistance is the ratio of physical transverse size to angular size
//   z : redshift
//   distance : [Mpc/rad]
func (cos WACDM) AngularDiameterDistance(z float64) (distance float64) {
	return cos.ComovingTransverseDistance(z) / (1 + z)
}

// ComovingTransverseDistance is the comoving distance at z as seen from z=0
//   z : redshift
//   distance : [Mpc/rad]
func (cos WACDM) ComovingTransverseDistance(z float64) (distance float64) {
	return cos.ComovingTransverseDistanceZ1Z2(0, z)
}

// ComovingTransverseDistanceZ1Z2 is the comoving distance at z2 as seen from z1
//   z : redshift
//   distance : [Mpc/rad]
func (cos WACDM) ComovingTransverseDistanceZ1Z2(z1, z2 float64) (distance float64) {
	comovingDistance := cos.ComovingDistanceZ1Z2(z1, z2)
	Ok0 := 1 - (cos.Om0 + cos.Ol0)
	if Ok0 == 0 {
		return comovingDistance
	}

	hubbleDistance := cos.HubbleDistance()
	return hubbleDistance / math.Sqrt(math.Abs(Ok0)) *
		math.Sinh(math.Sqrt(math.Abs(Ok0))*comovingDistance/hubbleDistance)
}

// HubbleDistance is the inverse of the Hubble parameter
//   distance : [Mpc]
func (cos WACDM) HubbleDistance() float64 {
	return SpeedOfLightKmS / cos.H0
}

// ComovingDistance is the distance that is constant with the Hubble flow
// expressed in the physical distance at z=0.
//   z: redshift
//   distance: [Mpc]
//
// As the scale factor a = 1/(1+z) increases from 0.5 to 1,
// two objects separated by a proper distance of 10 Mpc at a=0.5 (z=1)
// will be separated by a proper distance of 2*10 Mpc at a=1.0 (z=0).
// The comoving distance between these objects is 20 Mpc.
func (cos WACDM) ComovingDistance(z float64) (distance float64) {
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
func (cos WACDM) ComovingDistanceZ1Z2Elliptic(z1, z2 float64) (distance float64) {
	s := math.Pow((1-cos.Om0)/cos.Om0, 1./3)
	prefactor := (SpeedOfLightKmS / cos.H0) * (1 / math.Sqrt(s*cos.Om0))
	return prefactor * (tElliptic(s/(1+z1)) - tElliptic(s/(1+z2)))
}

// ComovingDistanceZ1Z2Integrate is the comoving distance between two z
// in a flat lambda CDM cosmology using fixed Gaussian quadrature integration.
//   z1 : redshift
//   z2 : redshift
//   distance : [Mpc]
func (cos WACDM) ComovingDistanceZ1Z2Integrate(z1, z2 float64) (distance float64) {
	n := 1000 // Integration will be n-point Gaussian quadrature
	return cos.HubbleDistance() * quad.Fixed(cos.Einv, z1, z2, n, nil, 0)
}

// ComovingDistanceZ1Z2 is the base function for calculation of comoving distances
// Here is where the choice of fundamental calculation method is made:
// Elliptic integral, quadrature integration, or analytic for special cases.
//   z1 : redshift
//   z2 : redshift
//   distance : [Mpc]
func (cos WACDM) ComovingDistanceZ1Z2(z1, z2 float64) (distance float64) {
	switch {
	// Test for Ol0==0 first so that (Om0, Ol0) = (1, 0)
	// is handled by the analytic solution
	// rather than the explicit integration.
	case cos.Ol0 == 0:
		return cos.ComovingDistanceOMZ1Z2(z1, z2)
	case (cos.WA == 0) && (cos.W0 == -1) && (cos.Om0 < 1):
		return cos.ComovingDistanceZ1Z2Elliptic(z1, z2)
	default:
		return cos.ComovingDistanceZ1Z2Integrate(z1, z2)
	}
}

// ComovingDistanceOM is the analytic case of Omega_total=Omega_M
//   z: redshift
//   distance: [Mpc]
func (cos WACDM) ComovingDistanceOM(z float64) (distance float64) {
	// Call the internal function that just takes direct arguments
	// with nothing passed via the struct.
	return comovingDistanceOM(z, cos.Om0, cos.H0)
}

// ComovingDistanceOMZ1Z2 is the analytic case of Omega_total=Omega_M
// for the distance between two redshifts.
//   z : redshift
//   time :[Gyr]
//
// This *Z1Z2 form exists to parallel the other versions
// and allow it to be a shortcut option in ComovingDistanceZ1Z2.
// Naively, it's twice as expensive to do this as (0, z2)
// But this is such a trivial calculation, it probably doesn't matter.
func (cos WACDM) ComovingDistanceOMZ1Z2(z1, z2 float64) (distance float64) {
	return comovingDistanceOM(z2, cos.Om0, cos.H0) -
		comovingDistanceOM(z1, cos.Om0, cos.H0)
}

// LookbackTime is the time from redshift 0 to z in Gyr.
//   z : redshift
//   time :[Gyr]
func (cos WACDM) LookbackTime(z float64) (time float64) {
	switch {
	case (cos.Ol0 == 0) && (0 < cos.Om0) && (cos.Om0 != 1):
		return cos.LookbackTimeOM(z)
	default:
		return cos.LookbackTimeIntegrate(z)
	}
}

// LookbackTimeIntegrate is the lookback time using explicit integration
//   z : redshift
//   time :[Gyr]
func (cos WACDM) LookbackTimeIntegrate(z float64) (time float64) {
	n := 1000 // Integration will be n-point Gaussian quadrature
	integrand := func(z float64) float64 { return cos.Einv(z) / (1 + z) }
	return hubbleTime(cos.H0) * quad.Fixed(integrand, 0, z, n, nil, 0)
}

// LookbackTimeOL is lookback time for dark-energy only Universe
//   z : redshift
//   time :[Gyr]
func (cos WACDM) LookbackTimeOL(z float64) (time float64) {
	return lookbackTimeOL(z, cos.Ol0, cos.H0)
}

// LookbackTimeOM is lookback time for matter only Universe
// All matter is non-relativistic.
//   z : redshift
//   time :[Gyr]
func (cos WACDM) LookbackTimeOM(z float64) (time float64) {
	return lookbackTimeOM(z, cos.Om0, cos.H0)
}

// Age is the time from redshift ∞ to z in Gyr.
//   z : redshift
//   time :[Gyr]
func (cos WACDM) Age(z float64) (time float64) {
	switch {
	case (cos.Ol0 == 0) && (0 < cos.Om0) && (cos.Om0 != 1):
		return cos.AgeOM(z)
	case (cos.W0 == -1.0) && (cos.Om0 == 0) && (0 < cos.Ol0) && (cos.Ol0 < 1):
		return cos.AgeOL(z)
	case (cos.W0 == -1.0) && (cos.Om0+cos.Ol0 == 1):
		return cos.AgeFlatLCDM(z)
	default:
		return cos.AgeIntegrate(z)
	}
}

// AgeFlatLCDM is the time from redshift ∞ to z
// in a flat LCDM cosmology.
//   z : redshift
//   time :[Gyr]
//
// Equation is in many sources.
// I took this from Thomas and Kantowski, 2000 PRD, 62, 103507.
func (cos WACDM) AgeFlatLCDM(z float64) (time float64) {
	return ageFlatLCDM(z, cos.Om0, cos.H0)
}

// AgeIntegrate is the time from redshift ∞ to z
// using explicit integration.
//   z : redshift
//   time :[Gyr]
//
// Basic integrand can be found in many texts
// I happened to copy this from
// Thomas and Kantowski, 2000, PRD, 62, 103507.  Eq. 1.
// Current implementation is fixed quadrature using mathext.integrate.quad.Fixed
func (cos WACDM) AgeIntegrate(z float64) (time float64) {
	n := 1000 // Integration will be n-point Gaussian quadrature
	integrand := func(z float64) float64 {
		denom := (1 + z) * cos.E(z)
		return 1 / denom
	}
	// When given math.Inf(), quad.Fixed automatically redefines variables
	// to successfully do the numerical integration.
	return hubbleTime(cos.H0) * quad.Fixed(integrand, z, math.Inf(1), n, nil, 0)
}

// AgeOL is the time from redshift ∞ to z
// with only constant dark energy and curvature.
//   z : redshift
//   time :[Gyr]
func (cos WACDM) AgeOL(z float64) (time float64) {
	return ageOL(z, cos.Ol0, cos.H0)
}

// AgeOL is the time from redshift ∞ to z
// with only non-relativistic matter and curvature.
//   z : redshift
//   time :[Gyr]
func (cos WACDM) AgeOM(z float64) (time float64) {
	return ageOM(z, cos.Om0, cos.H0)
}

// E is the Hubble parameter as a fraction of its present value.
//   z : redshift
//   ez : fractional Hubble parameter
// E.g., Hogg arXiv:9905116  Eq. 14
// Linder, 2003, PhRvL, 90, 130, Eq. 5, 7
func (cos WACDM) E(z float64) (ez float64) {
	oR := cos.Ogamma0 + cos.Onu0
	var deScale float64
	switch {
	case (cos.W0 == -1) && (cos.WA == 0):
		deScale = 1
	case cos.WA == 0:
		deScale = math.Pow(1+z, 3*(1+cos.W0))
	default:
		deScale = math.Pow(1+z, 3*(1+cos.W0+cos.WA)) * math.Exp(-3*cos.WA*z/(1+z))
	}
	Ok0 := 1 - (cos.Om0 + cos.Ol0)
	ez = math.Sqrt((1+z)*(1+z)*(1+z)*(1+z)*oR + (1+z)*(1+z)*(1+z)*cos.Om0 +
		(1+z)*(1+z)*Ok0 + cos.Ol0*deScale)
	return ez
}

// Einv is the inverse Hubble parameter
// Implementation is just to return E(z)
//   z : redshift
//   invEz : 1 / (fractional Hubble parameter)
func (cos WACDM) Einv(z float64) (invEz float64) {
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
