package cosmo

import (
	"fmt"
)

func ExampleWCDM_DistanceModulus() {
	var z_vec, obs_vec []float64
	cos := WCDM{Om0: 0.3, Ol0: 0.7, W0: -1.2, H0: 70, Tcmb0: 0.}

	z_vec = []float64{0.5, 1.0, 2.0, 3.0}
	obs_vec = make([]float64, len(z_vec))
	for i, z := range z_vec {
		obs_vec[i] = cos.DistanceModulus(z)
	}

	fmt.Println(obs_vec)
	// Output:
	// [42.26118542154089 44.10023765554372 45.95719725271018 47.026111928689645]
}

func ExampleWCDM_LuminosityDistance() {
	var z_vec, obs_vec []float64
	cos := WCDM{Om0: 0.3, Ol0: 0.7, W0: -1, H0: 70, Tcmb0: 0.}

	z_vec = []float64{0.5, 1.0, 2.0, 3.0}
	obs_vec = make([]float64, len(z_vec))
	for i, z := range z_vec {
		obs_vec[i] = cos.LuminosityDistance(z)
	}

	fmt.Println(obs_vec)
	// Output:
	// [2832.938093900105 6607.65761177494 15539.58622322811 25422.74174518986]
}

func ExampleWCDM_AngularDiameterDistance() {
	var z_vec, obs_vec []float64
	cos := WCDM{Om0: 0.3, Ol0: 0.7, W0: -1, H0: 70, Tcmb0: 0.}

	z_vec = []float64{0.5, 1.0, 2.0, 3.0}
	obs_vec = make([]float64, len(z_vec))
	for i, z := range z_vec {
		obs_vec[i] = cos.AngularDiameterDistance(z)
	}

	fmt.Println(obs_vec)
	// Output:
	// [1259.0835972889354 1651.914402943735 1726.6206914697902 1588.9213590743661]
}
