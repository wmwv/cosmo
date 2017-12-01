package cosmo

import (
	"reflect"
	"testing"
)

func benchmarkEN(n int, b *testing.B) {
	cos := Cosmology{Om0: 0.27, Ol0: 0.73, Ok0: 0., H0: 70, w0: -1.0, Tcmb0: 0.}

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

func BenchmarkEN(b *testing.B) {
	benchmarkEN(10000, b)
}

func BenchmarkENdistance(b *testing.B) {
	benchmarkNdistance(10000, "E", b)
}

func BenchmarkE(b *testing.B) {
	cos := Cosmology{Om0: 0.27, Ol0: 0.73, Ok0: 0., H0: 70, w0: -1.0, Tcmb0: 0.}
	z := 1.0
	for i := 0; i < b.N; i++ {
		cos.E(z)
	}
}

func BenchmarkEinv(b *testing.B) {
	cos := Cosmology{Om0: 0.27, Ol0: 0.73, Ok0: 0., H0: 70, w0: -1.0, Tcmb0: 0.}
	z := 1.0
	for i := 0; i < b.N; i++ {
		cos.Einv(z)
	}
}

// benchmarkDistanceOM is a helper function to be called by specific benchmarks
//   for an Omega_Lambda = 0 cosmology
func benchmarkDistanceOM(distFunc string, b *testing.B) {
    cos := Cosmology{Om0: 0.27, Ol0: 0., Ok0: 0., H0: 70, w0: -1.0, Tcmb0: 0.}
    z := 1.0

    funcToTest := reflect.ValueOf(&cos).MethodByName(distFunc)
    for i := 0; i < b.N; i++ {
        funcToTest.Call([]reflect.Value{reflect.ValueOf(z)})
    }
}

// benchmarkDistance is a helper function to be called by specific benchmarks
func benchmarkDistance(distFunc string, b *testing.B) {
	cos := Cosmology{Om0: 0.27, Ol0: 0.73, Ok0: 0., H0: 70, w0: -1.0, Tcmb0: 0.}
	z := 1.0

	funcToTest := reflect.ValueOf(&cos).MethodByName(distFunc)
	for i := 0; i < b.N; i++ {
		funcToTest.Call([]reflect.Value{reflect.ValueOf(z)})
	}
}

// benchmarkNdistance is a helper function to be called by specific benchmarks
func benchmarkNdistance(n int, distFunc string, b *testing.B) {
	cos := Cosmology{Om0: 0.27, Ol0: 0.73, Ok0: 0., H0: 70, w0: -1.0, Tcmb0: 0.}
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

func BenchmarkComovingDistance(b *testing.B) {
	benchmarkDistance("ComovingDistance", b)
}

func BenchmarkComovingTransverseDistance(b *testing.B) {
	benchmarkDistance("ComovingTransverseDistance", b)
}

func BenchmarkLuminosityDistance(b *testing.B) {
	benchmarkDistance("LuminosityDistance", b)
}

func BenchmarkNComovingDistance(b *testing.B) {
	benchmarkNdistance(10000, "ComovingDistance", b)
}

func BenchmarkNLuminosityDistance(b *testing.B) {
	benchmarkNdistance(10000, "LuminosityDistance", b)
}

func BenchmarkNE(b *testing.B) {
	benchmarkNdistance(10000, "E", b)
}

func BenchmarkComovingDistanceOM(b *testing.B) {
    benchmarkDistanceOM("ComovingDistance", b)
}
