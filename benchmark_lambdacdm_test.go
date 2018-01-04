package cosmo

import (
	"reflect"
	"testing"
)

func benchmarkLambdaCDMEN(n int, b *testing.B) {
	cos := LambdaCDM{Om0: 0.27, Ol0: 0.73, H0: 70}

	var z float64
	z_max := 1.0
	step := z_max / float64(n)
	for i := 0; i < b.N; i++ {
		for j := 0; j < n; j++ {
			z = 0.001 + step*float64(j)
			cos.E(z)
		}
	}
}

func BenchmarkLambdaCDMEN(b *testing.B) {
	benchmarkLambdaCDMEN(10000, b)
}

func BenchmarkLambdaCDMENdistance(b *testing.B) {
	benchmarkLambdaCDMNdistance(10000, "E", b)
}

func BenchmarkLambdaCDME(b *testing.B) {
	cos := LambdaCDM{Om0: 0.27, Ol0: 0.73, H0: 70}
	z := 1.0
	for i := 0; i < b.N; i++ {
		cos.E(z)
	}
}

func BenchmarkLambdaCDMEinv(b *testing.B) {
	cos := LambdaCDM{Om0: 0.27, Ol0: 0.73, H0: 70}
	z := 1.0
	for i := 0; i < b.N; i++ {
		cos.Einv(z)
	}
}

// benchmarkLambdaCDMDistanceOM is a helper function to be called by specific benchmarkLambdaCDMs
//   for an Omega_Lambda = 0 cosmology
func benchmarkLambdaCDMDistanceOM(distFunc string, b *testing.B) {
	cos := LambdaCDM{Om0: 0.27, Ol0: 0., H0: 70}
	z := 1.0

	funcToTest := reflect.ValueOf(&cos).MethodByName(distFunc)
	for i := 0; i < b.N; i++ {
		funcToTest.Call([]reflect.Value{reflect.ValueOf(z)})
	}
}

// benchmarkLambdaCDMDistance is a helper function to be called by specific benchmarkLambdaCDMs
func benchmarkLambdaCDMDistance(distFunc string, b *testing.B) {
	cos := LambdaCDM{Om0: 0.27, Ol0: 0.73, H0: 70}
	z := 1.0

	funcToTest := reflect.ValueOf(&cos).MethodByName(distFunc)
	for i := 0; i < b.N; i++ {
		funcToTest.Call([]reflect.Value{reflect.ValueOf(z)})
	}
}

// benchmarkLambdaCDMNdistance is a helper function to be called by specific benchmarkLambdaCDMs
func benchmarkLambdaCDMNdistance(n int, distFunc string, b *testing.B) {
	cos := LambdaCDM{Om0: 0.27, Ol0: 0.73, H0: 70}
	funcToTest := reflect.ValueOf(&cos).MethodByName(distFunc)
	var z float64
	z_max := 1.0
	step := z_max / float64(n)
	for i := 0; i < b.N; i++ {
		for j := 0; j < n; j++ {
			z = 0.001 + step*float64(j)
			funcToTest.Call([]reflect.Value{reflect.ValueOf(z)})
		}
	}
}

func BenchmarkLambdaCDMComovingDistance(b *testing.B) {
	benchmarkLambdaCDMDistance("ComovingDistance", b)
}

func BenchmarkLambdaCDMComovingTransverseDistance(b *testing.B) {
	benchmarkLambdaCDMDistance("ComovingTransverseDistance", b)
}

func BenchmarkLambdaCDMLuminosityDistance(b *testing.B) {
	benchmarkLambdaCDMDistance("LuminosityDistance", b)
}

func BenchmarkLambdaCDMLookbackTime(b *testing.B) {
	benchmarkLambdaCDMDistance("LookbackTime", b)
}

func BenchmarkLambdaCDMNComovingDistance(b *testing.B) {
	benchmarkLambdaCDMNdistance(10000, "ComovingDistance", b)
}

func BenchmarkLambdaCDMNLuminosityDistance(b *testing.B) {
	benchmarkLambdaCDMNdistance(10000, "LuminosityDistance", b)
}

func BenchmarkLambdaCDMNE(b *testing.B) {
	benchmarkLambdaCDMNdistance(10000, "E", b)
}

func BenchmarkLambdaCDMComovingDistanceOM(b *testing.B) {
	benchmarkLambdaCDMDistanceOM("ComovingDistance", b)
}

func BenchmarkLambdaCDMLookbackTimeOM(b *testing.B) {
	benchmarkLambdaCDMDistanceOM("LookbackTime", b)
}
