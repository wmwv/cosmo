// Package cosmogo implements basic cosmology calculations in Go
// Provides E, LuminosityDistance, DistanceModulus

// Based on code in astropy.cosmology
//   http://docs.astropy.org/en/stable/_modules/astropy/cosmology
// And David Hogg's
//   https://arxiv.org/abs/astro-ph/9905116
package cosmo

import (
	"gonum.org/v1/gonum/integrate/quad"
	"gonum.org/v1/gonum/mathext"
	"math"
)

const SpeedOfLightKmS = 299792.458 // km/s

// Cosmology stores the key information needed for a given cosmology
type Cosmology struct {
	Om0     float64 // Matter Density at z=0
	Ol0     float64 // Vacuum Energy density Lambda at z=0
	Ok0     float64 // Curvature Density at z=0
	H0      float64 // Hubble constant at z=0.  [km/s/Mpc]
	w0      float64 // Dark energy equation-of-state parameter
	Ogamma0 float64 // Photon density
	Onu0    float64 // Neutrino density
	Tcmb0   float64 // Temperature of the CMB at z=0.  [K]
	//    nuToPhotonDensity float64 // Neutrino density / photon density
}

// DistanceModulus returns the magnitude difference between 1 Mpc and
// the luminosity distance for the given z.
func (cos *Cosmology) DistanceModulus(z float64) (distance float64) {
	return 5*math.Log10(cos.LuminosityDistance(z)) + 25
}

func (cos *Cosmology) LuminosityDistance(z float64) (distance float64) {
	return (1 + z) * cos.ComovingTransverseDistance(z)
}

// ComovingTransversedistance is the comoving distance at z as seen from z=0
// Returns scalar in Megaparsecs
func (cos *Cosmology) ComovingTransverseDistance(z float64) (distance float64) {
	return cos.ComovingTransverseDistanceZ1Z2(0, z)
}

// ComovingTransverseDistanceZ1Z2 is the comoving distance at z2 as seen from z1
// Returns scalar in Megaparsecs
func (cos *Cosmology) ComovingTransverseDistanceZ1Z2(z1, z2 float64) (distance float64) {
	comovingDistance := cos.ComovingDistanceZ1Z2(z1, z2)
	if cos.Ok0 == 0 {
		return comovingDistance
	}

	hubbleDistance := cos.HubbleDistance()
	return hubbleDistance /
		math.Sinh(math.Sqrt(cos.Ok0)*comovingDistance/hubbleDistance)
}

func (cos *Cosmology) HubbleDistance() float64 {
	return SpeedOfLightKmS / cos.H0
}

func (cos *Cosmology) ComovingDistance(z float64) (distance float64) {
	return cos.ComovingDistanceZ1Z2(0, z)
}

// ComovingDistanceElliptic calculates the comoving distance from 0 to z
//   in a flat lambda CDM cosmology using elliptic integrals.
func (cos *Cosmology) ComovingDistanceElliptic(z float64) (distance float64) {
	return cos.ComovingDistanceZ1Z2Elliptic(0, z)
}

// ComovingDistanceElliptic calculates the comoving distance between two z
//   in a flat lambda CDM cosmology using elliptic integrals.
func (cos *Cosmology) ComovingDistanceZ1Z2Elliptic(z1, z2 float64) (distance float64) {
	s := math.Pow((1-cos.Om0)/cos.Om0, 1./3)
	prefactor := (SpeedOfLightKmS / cos.H0) * (1 / math.Sqrt(s*cos.Om0))
	return prefactor * (TElliptic(s/(1+z1)) - TElliptic(s/(1+z2)))
}

// TElliptic uses elliptic integral of the first kind
//   basic integral at the heart of cosmological distance calculations
// gonum.org/v1/mathext/EllipticF (Legendre form) reformats
//   and calls EllipticRF (Carlson form)
func TElliptic(x float64) float64 {
	numer := 1 + (1-math.Sqrt(3))*x
	denom := 1 + (1+math.Sqrt(3))*x
	phi := math.Acos(numer / denom)
	// math.Sqrt is ~35 times faster than math.Cos
	// Although I would expect these fixed constants
	//  to be calculated at compile time, I don't know if that happens.
	// k := math.Cos(math.Pi / 12)
	m := 0.5 + math.Sqrt(3)/4 // = cos(pi/12)**2
	return math.Pow(3, -1./4) * mathext.EllipticF(phi, m)
}

func (cos *Cosmology) ComovingDistanceZ1Z2Integrate(z1, z2 float64) (distance float64) {
	n := 1000 // Integration will be n-point Gaussian quadrature
	return cos.HubbleDistance() * quad.Fixed(cos.Einv, z1, z2, n, nil, 0)
}

func (cos *Cosmology) ComovingDistanceZ1Z2(z1, z2 float64) (distance float64) {
	switch {
	case cos.Om0 < 1:
		return cos.ComovingDistanceZ1Z2Elliptic(z1, z2)
	default:
		return cos.ComovingDistanceZ1Z2Integrate(z1, z2)
	}
}

// E calculates the Hubble parameter as a fraction of its present value
// E.g., Hogg arXiv:9905116  Eq. 14
//   E(z) = int_0^z (...) dz
// Integration is done via gonum.quad
func (cos *Cosmology) E(z float64) (ez float64) {
	oR := cos.Ogamma0 + cos.Onu0
	var deScale float64
	// TODO
	// Consider an if or switch on the value of cos.w0
	// Do performance testing to see in what circumstances it matters.
	switch cos.w0 {
	case -1:
		deScale = 1
	default:
		deScale = math.Pow(1+z, 3*(1+cos.w0))
	}
	ez = math.Sqrt((1+z)*(1+z)*((oR*(1+z)+cos.Om0)*(1+z)+cos.Ok0) + cos.Ol0*deScale)
	return ez
}

// inv_E is the thing generally integrated for calculating distances
// and related quantities
func (cos *Cosmology) Einv(z float64) (invEz float64) {
	// TODO: Is 1/Sqrt() notably slower than Pow(-0.5)?
	// No.
	// Pow(-0.5) is in fact implemented as 1/Sqrt() in math.pow.go
	// func pow(x, y float64) float64 {
	//    [...]
	// case y == -0.5:
	//    return 1 / Sqrt(x)
	return 1 / cos.E(z)
}

// Evec returns vectorized form of 'E'.
// I haven't figured out whether this is useful or makes sense
// in a Go framework, which looping over functions is more expected.
// thank in IDL, Matlab, or Python numpy+scipy worlds.
func (cos *Cosmology) Evec(z []float64) (ez []float64) {
	for _, z := range z {
		ez = append(ez, cos.E(z))
	}
	return ez
}
