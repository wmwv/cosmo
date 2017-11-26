// Package cosmogo implements basic cosmology calculations in Go
// Provides E, LuminosityDistance, DistanceModulus

// Based on code in astropy.cosmology
//   http://docs.astropy.org/en/stable/_modules/astropy/cosmology
// And David Hogg's
//   https://arxiv.org/abs/astro-ph/9905116
package cosmo

import (
	"gonum.org/v1/gonum/integrate/quad"
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

func (cos *Cosmology) ComovingDistanceZ1Z2(z1, z2 float64) (distance float64) {
	n := 10000 // integration points
	return cos.HubbleDistance() * quad.Fixed(cos.Einv, z1, z2, n, nil, 0)
}

// E calculates the Hubble parameter as a fraction of its present value
// E.g., Hogg arXiv:9905116  Eq. 14
//   E(z) = int_0^z (...) dz
// Integration is done via gonum.quad
func (cos *Cosmology) E(z float64) (ez float64) {
	oR := cos.Ogamma0 + cos.Onu0
	// TODO
	// Consider an if or switch on the value of cos.w0
	// Do performance testing to see in what circumstances it matters.
	ez = math.Sqrt((1+z)*(1+z)*
		((oR*(1+z)+cos.Om0)*(1+z)+cos.Ok0) +
		cos.Ol0*math.Pow(1+z, 3*(1+cos.w0)))
	return ez
}

// inv_E is the thing generally integrated for calculating distances
// and related quantities
func (cos *Cosmology) Einv(z float64) (invEz float64) {
	// TODO: Is 1/Sqrt() notable slower than Pow(-0.5)?
	// If so, then it's worth writing out inv_E again here.
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
