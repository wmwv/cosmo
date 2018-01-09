package cosmo

import (
	"reflect"
	"testing"
)

func benchmarkWACDMEN(n int, b *testing.B) {
	cos := WACDM{H0: 70, Om0: 0.2, Ol0: 0.7, W0: -1.2, WA: 2}

	var z float64
	zMax := 1.0
	step := zMax / float64(n)
	for i := 0; i < b.N; i++ {
		for j := 0; j < n; j++ {
			z = 0.001 + step*float64(j)
			cos.E(z)
		}
	}
}

func BenchmarkWACDMEN(b *testing.B) {
	benchmarkWACDMEN(10000, b)
}

func BenchmarkWACDMENdistance(b *testing.B) {
	benchmarkWACDMNdistance(10000, "E", b)
}

func BenchmarkWACDME(b *testing.B) {
	cos := WACDM{H0: 70, Om0: 0.2, Ol0: 0.7, W0: -1.2, WA: 2}
	z := 1.0
	for i := 0; i < b.N; i++ {
		cos.E(z)
	}
}

func BenchmarkWACDMEinv(b *testing.B) {
	cos := WACDM{H0: 70, Om0: 0.2, Ol0: 0.7, W0: -1.2, WA: 2}
	z := 1.0
	for i := 0; i < b.N; i++ {
		cos.Einv(z)
	}
}

// benchmarkWACDMDistanceOM is a helper function to be called by specific benchmarkWACDMs
//   for an Omega_Lambda = 0 cosmology
func benchmarkWACDMDistanceOM(distFunc string, b *testing.B) {
	cos := WACDM{H0: 70, Om0: 0.3, Ol0: 0, W0: -1.2, WA: 2}
	z := 1.0

	funcToTest := reflect.ValueOf(&cos).MethodByName(distFunc)
	for i := 0; i < b.N; i++ {
		funcToTest.Call([]reflect.Value{reflect.ValueOf(z)})
	}
}

// benchmarkWACDMDistance is a helper function to be called by specific benchmarkWACDMs
func benchmarkWACDMDistance(distFunc string, b *testing.B) {
	cos := WACDM{H0: 70, Om0: 0.2, Ol0: 0.7, W0: -1.2}
	z := 1.0

	funcToTest := reflect.ValueOf(&cos).MethodByName(distFunc)
	for i := 0; i < b.N; i++ {
		funcToTest.Call([]reflect.Value{reflect.ValueOf(z)})
	}
}

// benchmarkWACDMNdistance is a helper function to be called by specific benchmarkWACDMs
func benchmarkWACDMNdistance(n int, distFunc string, b *testing.B) {
	cos := WACDM{H0: 70, Om0: 0.2, Ol0: 0.7, W0: -1.2}
	funcToTest := reflect.ValueOf(&cos).MethodByName(distFunc)
	var z float64
	zMax := 1.0
	step := zMax / float64(n)
	for i := 0; i < b.N; i++ {
		for j := 0; j < n; j++ {
			z = 0.001 + step*float64(j)
			funcToTest.Call([]reflect.Value{reflect.ValueOf(z)})
		}
	}
}

func BenchmarkWACDMComovingDistance(b *testing.B) {
	benchmarkWACDMDistance("ComovingDistance", b)
}

func BenchmarkWACDMComovingTransverseDistance(b *testing.B) {
	benchmarkWACDMDistance("ComovingTransverseDistance", b)
}

func BenchmarkWACDMLuminosityDistance(b *testing.B) {
	benchmarkWACDMDistance("LuminosityDistance", b)
}

func BenchmarkWACDMLookbackTime(b *testing.B) {
	benchmarkWACDMDistance("LookbackTime", b)
}

func BenchmarkWACDMNComovingDistance(b *testing.B) {
	benchmarkWACDMNdistance(10000, "ComovingDistance", b)
}

func BenchmarkWACDMNLuminosityDistance(b *testing.B) {
	benchmarkWACDMNdistance(10000, "LuminosityDistance", b)
}

func BenchmarkWACDMNE(b *testing.B) {
	benchmarkWACDMNdistance(10000, "E", b)
}

func BenchmarkWACDMComovingDistanceOM(b *testing.B) {
	benchmarkWACDMDistanceOM("ComovingDistance", b)
}

func BenchmarkWACDMLookbackTimeOM(b *testing.B) {
	benchmarkWACDMDistanceOM("LookbackTime", b)
}
