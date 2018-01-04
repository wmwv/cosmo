package cosmo

import (
	"gonum.org/v1/gonum/integrate/quad"
	"math"
)

// FlatLCDM provides cosmological distances, age, and look-back time
// for a flat Universe with a cosmological constant:
// matter, dark energy, and no curvature,
// with a w=-1 equation-of-state parameter for dark energy
//
// No radiation pressure or neutrino contributions.
type FlatLCDM struct {
	Om0     float64 // Matter Density at z=0
	H0      float64 // Hubble constant at z=0.  [km/s/Mpc]
	W0      float64 // Dark energy equation-of-state parameter
	Ogamma0 float64 // Photon density
	Onu0    float64 // Neutrino density
}

// Tcmb0   float64 // Temperature of the CMB at z=0.  [K]
// nuToPhotonDensity float64 // Neutrino density / photon density

// DistanceModulus is the magnitude difference between 1 Mpc and
// the luminosity distance for the given z.
func (cos FlatLCDM) DistanceModulus(z float64) (distanceModulusMag float64) {
	return 5*math.Log10(cos.LuminosityDistance(z)) + 25
}

// LuminosityDistance is the radius of effective sphere over which the light has spread out
func (cos FlatLCDM) LuminosityDistance(z float64) (distanceMpc float64) {
	return (1 + z) * cos.ComovingTransverseDistance(z)
}

// AngularDistance is the ratio of physical transverse size to angular size
func (cos FlatLCDM) AngularDiameterDistance(z float64) (distanceMpcRad float64) {
	return cos.ComovingTransverseDistance(z) / (1 + z)
}

// ComovingTransverseDistance is the comoving distance at z as seen from z=0
func (cos FlatLCDM) ComovingTransverseDistance(z float64) (distanceMpcRad float64) {
	return cos.ComovingTransverseDistanceZ1Z2(0, z)
}

// ComovingTransverseDistanceZ1Z2 is the comoving distance at z2 as seen from z1
func (cos FlatLCDM) ComovingTransverseDistanceZ1Z2(z1, z2 float64) (distanceMpcRad float64) {
	comovingDistance := cos.ComovingDistanceZ1Z2(z1, z2)
	return comovingDistance
}

// HubbleDistance is the inverse of the Hubble parameter times the speed of light.
func (cos FlatLCDM) HubbleDistance() (distanceMpc float64) {
	return SpeedOfLightKmS / cos.H0
}

// ComovingDistance is the distance that is constant with the Hubble flow
// expressed in the physical distance at z=0.
//
// As the scale factor a = 1/(1+z) increases from 0.5 to 1,
// two objects separated by a proper distance of 10 Mpc at a=0.5 (z=1)
// will be separated by a proper distance of 10*2 Mpc at a=1.0 (z=0)
// The comoving distance between these objects is 20 Mpc
func (cos FlatLCDM) ComovingDistance(z float64) (distanceMpc float64) {
	return cos.ComovingDistanceZ1Z2(0, z)
}

// ComovingDistanceZ1Z2Elliptic is the comoving distance between two redshifts
// in a flat lambda CDM cosmology using elliptic integrals.
//
// See
//     Feige, 1992, Astron. Nachr., 313, 139.
//     Eisenstein, 1997, arXiv:9709054v2
//     Mészáros & Řípai 2013, A&A, 556, A13.
// and a useful summary in
//     Baes, Camps, Van De Putte, 2017, MNRAS, 468, 927.
func (cos FlatLCDM) comovingDistanceZ1Z2Elliptic(z1, z2 float64) (distanceMpc float64) {
	s := math.Pow((1-cos.Om0)/cos.Om0, 1./3)
	prefactor := (SpeedOfLightKmS / cos.H0) * (1 / math.Sqrt(s*cos.Om0))
	return prefactor * (tElliptic(s/(1+z1)) - tElliptic(s/(1+z2)))
}

// ComovingDistanceZ1Z2Integrate is the comoving distance between two z
// in a flat lambda CDM cosmology using fixed Gaussian quadrature integration.
func (cos FlatLCDM) comovingDistanceZ1Z2Integrate(z1, z2 float64) (distanceMpc float64) {
	n := 1000 // Integration will be n-point Gaussian quadrature
	return cos.HubbleDistance() * quad.Fixed(cos.Einv, z1, z2, n, nil, 0)
}

// ComovingDistanceZ1Z2 is the base function for calculation of comoving distances
//
// Here is where the choice of fundamental calculation method is made:
// Elliptic integral, quadrature integration, or analytic for special cases.
func (cos FlatLCDM) ComovingDistanceZ1Z2(z1, z2 float64) (distanceMpc float64) {
	switch {
	case cos.Om0 < 1:
		return cos.comovingDistanceZ1Z2Elliptic(z1, z2)
	default:
		return cos.comovingDistanceZ1Z2Integrate(z1, z2)
	}
}

// LookbackTime is the time from redshift 0 to z.
func (cos FlatLCDM) LookbackTime(z float64) (timeGyr float64) {
	return cos.lookbackTimeIntegrate(z)
}

// LookbackTimeIntegrate is the look-back time using explicit integration
func (cos FlatLCDM) lookbackTimeIntegrate(z float64) (timeGyr float64) {
	n := 1000 // Integration will be n-point Gaussian quadrature
	integrand := func(z float64) float64 { return cos.Einv(z) / (1 + z) }
	return hubbleTime(cos.H0) * quad.Fixed(integrand, 0, z, n, nil, 0)
}

// AgeFlatLCDM is the time from redshift ∞ to z
// in a flat LCDM cosmology.
func (cos FlatLCDM) Age(z float64) (timeGyr float64) {
	return ageFlatLCDM(z, cos.Om0, cos.H0)
}

// ageIntegrate is the time from redshift ∞ to z
// using explicit integration.
//
// Age is analytic in a Flat LCDM Universe.
// This function exists for consistency testing.
// The basic integrand can be found in many texts.
// I happened to copy this from
// Thomas and Kantowski, 2000, PRD, 62, 103507.  Eq. 1.
// Current implementation is fixed quadrature using mathext.integrate.quad.Fixed
func (cos FlatLCDM) ageIntegrate(z float64) (timeGyr float64) {
	n := 1000 // Integration will be n-point Gaussian quadrature
	Ol0 := 1 - cos.Om0
	integrand := func(z float64) float64 {
		denom := (1 + z) * math.Sqrt((1+z)*(1+z)*(1+cos.Om0*z)-z*(2+z)*Ol0)
		return 1 / denom
	}
	// When given math.Inf(), quad.Fixed automatically redefines variables
	// to successfully do the numerical integration.
	return hubbleTime(cos.H0) * quad.Fixed(integrand, z, math.Inf(1), n, nil, 0)
}

// E is the Hubble parameter as a fraction of its present value.
// E.g., Hogg arXiv:9905116  Eq. 14
func (cos FlatLCDM) E(z float64) (fractionalHubbleParameter float64) {
	oR := cos.Ogamma0 + cos.Onu0
	deScale := 1.0
	Ol0 := 1 - cos.Om0
	// I find it easier to explicitly // see the powers of 4, 3, and 0
	// The compiler is good enough that this operation or
	// 1+z vs. opz doesn't matter
	return math.Sqrt((1+z)*(1+z)*(1+z)*(1+z)*oR + (1+z)*(1+z)*(1+z)*cos.Om0 +
		Ol0*deScale)
}

// Einv is the inverse Hubble parameter
func (cos FlatLCDM) Einv(z float64) (invFractionalHubbleParameter float64) {
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
