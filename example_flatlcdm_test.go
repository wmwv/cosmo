package cosmo

import (
	"fmt"
)

// Calculated via
//   from astropy.cosmology import FlatLCDM
//   z = np.asarray([0.5, 1.0, 2.0, 3.0])
//   FlatLCDM(70, 0.3, 0.7).distmod(z)
//   FlatLCDM(70, 0.3, 0.7).luminosity_distance(z)
//   FlatLCDM(70, 0.3, 0.7).angular_diameter_distance(z)

func ExampleFlatLCDM() {
	cos := FlatLCDM{H0: 70, Om0: 0.3}

	zVec := []float64{0.5, 1.0, 2.0, 3.0}
	distmodVec := make([]float64, len(zVec))
	lumdistVec := make([]float64, len(zVec))
	angdistVec := make([]float64, len(zVec))
	for i, z := range zVec {
		distmodVec[i] = cos.DistanceModulus(z)
		lumdistVec[i] = cos.LuminosityDistance(z)
		angdistVec[i] = cos.AngularDiameterDistance(z)
	}

	fmt.Println(cos)
	fmt.Println("Ok0: ", cos.Ok0())
	fmt.Println("Distance Modulus [mag]")
	fmt.Println(distmodVec)
	fmt.Println("Luminosity Distance [Mpc]")
	fmt.Println(lumdistVec)
	fmt.Println("Angular Diameter Distance [Mpc]")
	fmt.Println(angdistVec)
	// Output:
	// FlatLCDM{H0: 70, Om0: 0.3}
	// Ok0:  0
	// Distance Modulus [mag]
	// [42.26118542154089 44.10023765554372 45.95719725271018 47.026111928689645]
	// Luminosity Distance [Mpc]
	// [2832.938093900105 6607.65761177494 15539.58622322811 25422.74174518986]
	// Angular Diameter Distance [Mpc]
	// [1259.0835972889354 1651.914402943735 1726.6206914697902 1588.9213590743661]
}
