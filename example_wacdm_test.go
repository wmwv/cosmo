package cosmo

import (
	"fmt"
)

// Calculated via
//   from astropy.cosmology import w0waCDM
//   z = np.asarray([0.5, 1.0, 2.0, 3.0])
//   w0waCDM(70, 0.3, 0.7, -0.8, 2.5).distmod(z)
//   w0waCDM(70, 0.3, 0.7, -0.8, 2.5).luminosity_distance(z)
//   w0waCDM(70, 0.3, 0.7, -0.8, 2.5).angular_diameter_distance(z)

func ExampleWACDM() {
	cos := WACDM{H0: 70, Om0: 0.3, Ol0: 0.7, W0: -0.8, WA: 2.5}

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
	// Output:
	// Distance Modulus [mag]
	// [42.07480332804884 43.731011211176536 45.31078970620773 46.17487505099648]
	// Luminosity Distance [Mpc]
	// [2599.9240753482904 5574.452795915061 11538.72814084889 17178.09539590076]
	// Angular Diameter Distance [Mpc]
	// [1155.521811265907 1393.6131989787652 1282.0809045387657 1073.6309622437975]
}
