package cosmo

import (
	"fmt"
)

// Calculated via
//   from astropy.cosmology import w0waCDM
//   z = np.asarray([0.5, 1.0, 2.0, 3.0])
//   wCDM(70, 0.3, 0.7, -1.2).distmod(z)
//   wCDM(70, 0.3, 0.7, -1.2).luminosity_distance(z)
//   wCDM(70, 0.3, 0.7, -1.2).angular_diameter_distance(z)

func ExampleWCDM() {
	cos := WCDM{H0: 70, Om0: 0.3, Ol0: 0.7, W0: -1.2}

	z_vec := []float64{0.5, 1.0, 2.0, 3.0}
	distmod_vec := make([]float64, len(z_vec))
	lumdist_vec := make([]float64, len(z_vec))
	angdist_vec := make([]float64, len(z_vec))
	for i, z := range z_vec {
		distmod_vec[i] = cos.DistanceModulus(z)
		lumdist_vec[i] = cos.LuminosityDistance(z)
		angdist_vec[i] = cos.AngularDiameterDistance(z)
	}

	fmt.Println(cos)
	fmt.Println("Distance Modulus [mag]")
	fmt.Println(distmod_vec)
	fmt.Println("Luminosity Distance [Mpc]")
	fmt.Println(lumdist_vec)
	fmt.Println("Angular Diameter Distance [Mpc]")
	fmt.Println(angdist_vec)
	// Output:
	// WCDM{H0: 70, Om0: 0.3, Ol0: 0.7, W0: -1.2}
	// Distance Modulus [mag]
	// [42.32710910996119 44.17957200628159 46.03118143998202 47.092287353314816]
	// Luminosity Distance [Mpc]
	// [2920.2620320966266 6853.531311400255 16078.157845948543 26209.423639506458]
	// Angular Diameter Distance [Mpc]
	// [1297.8942364873897 1713.3828278500637 1786.4619828831712 1638.0889774691536]
}
