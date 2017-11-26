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

// Cosmology stores the key information needed for a given cosmology
type Cosmology struct {
	Om0     float64 // Matter Density at z=0
	Ol0     float64 // Vacuum Energy density Lambda at z=0
	Ok0     float64 // Curvature Density at z=0
	H0      float64 // Hubble constant at z=0.  km/s/Mpc
	w0      float64 // Dark energy equation-of-state parameter
	Ogamma0 float64 // Photon density
	Onu0    float64 // Neutrino density
	Tcmb0   float64 // Temperature of the CMb at z=0.
	//    nuToPhotonDensity float64 // Neutrino density / photon density
}

//   E(z) = int_0^z (...) dz
// Integration is done via gonum.quad
func (cos *Cosmology) somethingelse(z float64) (Ez float64) {
	n := 10000 // integration points
	Ez = quad.Fixed(cos.E, 0, z, n, nil, 0)
	return Ez
}

// E calculates the Hubble parameter as a fraction of its present value
// E.g., Hogg arXiv:9905116  Eq. 14
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
