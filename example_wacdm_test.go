package cosmo

import (
	"fmt"
)

func ExampleWACDM() {
	cos := WACDM{Om0: 0.3, Ol0: 0.7, W0: -0.8, WA: 2.5, H0: 70, Tcmb0: 0.}

	z_vec := []float64{0.5, 1.0, 2.0, 3.0}
	distmod_vec := make([]float64, len(z_vec))
	lumdist_vec := make([]float64, len(z_vec))
	angdist_vec := make([]float64, len(z_vec))
	for i, z := range z_vec {
		distmod_vec[i] = cos.DistanceModulus(z)
		lumdist_vec[i] = cos.LuminosityDistance(z)
		angdist_vec[i] = cos.AngularDiameterDistance(z)
	}

	fmt.Println("Distance Modulus [mag]")
	fmt.Println(distmod_vec)
	fmt.Println("Luminosity Distance [Mpc]")
	fmt.Println(lumdist_vec)
	fmt.Println("Angular Diameter Distance [Mpc]")
	fmt.Println(angdist_vec)
	// Calculated via
	//   from astropy.cosmology import w0waCDM
	//   z = np.asarray([0.5, 1.0, 2.0, 3.0])
	//   w0waCDM(70, 0.3, 0.7, -0.8, 2.5).distmod(z)
	//   w0waCDM(70, 0.3, 0.7, -0.8, 2.5).luminosity_distance(z)
	//   w0waCDM(70, 0.3, 0.7, -0.8, 2.5).angular_diameter_distance(z)
	// Output:
	// Distance Modulus [mag]
	// [42.07480333, 43.73101121, 45.31078971, 46.17487505]
	// Luminosity Distance [Mpc]
	// [2599.92407535,  5574.45279592, 11538.72814085, 17178.0953959]
	// Angular Diameter Distance [Mpc]
	// [1155.52181127, 1393.61319898, 1282.08090454, 1073.63096224]
}
