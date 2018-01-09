package cosmo

import (
	"reflect"
	"testing"
)

func benchmarkWCDMEN(n int, b *testing.B) {
	cos := WCDM{H0: 70, Om0: 0.2, Ol0: 0.7, W0: -1.2}

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

func BenchmarkWCDMEN(b *testing.B) {
	benchmarkWCDMEN(10000, b)
}

func BenchmarkWCDMENdistance(b *testing.B) {
	benchmarkWCDMNdistance(10000, "E", b)
}

func BenchmarkWCDME(b *testing.B) {
	cos := WCDM{H0: 70, Om0: 0.2, Ol0: 0.7, W0: -1.2}
	z := 1.0
	for i := 0; i < b.N; i++ {
		cos.E(z)
	}
}

func BenchmarkWCDMEinv(b *testing.B) {
	cos := WCDM{H0: 70, Om0: 0.2, Ol0: 0.7, W0: -1.2}
	z := 1.0
	for i := 0; i < b.N; i++ {
		cos.Einv(z)
	}
}

// benchmarkWCDMDistanceFlat is a helper function to be called by specific benchmarkWCDMs
//   for an Omega_Lambda = 0 cosmology
func benchmarkWCDMDistanceFlat(distFunc string, b *testing.B) {
	cos := WCDM{H0: 70, Om0: 0.3, Ol0: 0.7, W0: -1}
	z := 1.0

	funcToTest := reflect.ValueOf(&cos).MethodByName(distFunc)
	for i := 0; i < b.N; i++ {
		funcToTest.Call([]reflect.Value{reflect.ValueOf(z)})
	}
}

// benchmarkWCDMDistanceLCDM is a helper function to be called by specific benchmarkWCDMs
//   for an Omega_Lambda = 0 cosmology
func benchmarkWCDMDistanceLCDM(distFunc string, b *testing.B) {
	cos := WCDM{H0: 70, Om0: 0.3, Ol0: 0.9, W0: -1}
	z := 1.0

	funcToTest := reflect.ValueOf(&cos).MethodByName(distFunc)
	for i := 0; i < b.N; i++ {
		funcToTest.Call([]reflect.Value{reflect.ValueOf(z)})
	}
}

// benchmarkWCDMDistanceOM is a helper function to be called by specific benchmarkWCDMs
//   for an Omega_Lambda = 0 cosmology
func benchmarkWCDMDistanceOM(distFunc string, b *testing.B) {
	cos := WCDM{H0: 70, Om0: 0.3, Ol0: 0, W0: -1.2}
	z := 1.0

	funcToTest := reflect.ValueOf(&cos).MethodByName(distFunc)
	for i := 0; i < b.N; i++ {
		funcToTest.Call([]reflect.Value{reflect.ValueOf(z)})
	}
}

// benchmarkWCDMDistance is a helper function to be called by specific benchmarkWCDMs
func benchmarkWCDMDistance(distFunc string, b *testing.B) {
	cos := WCDM{H0: 70, Om0: 0.2, Ol0: 0.7, W0: -1.2}
	z := 1.0

	funcToTest := reflect.ValueOf(&cos).MethodByName(distFunc)
	for i := 0; i < b.N; i++ {
		funcToTest.Call([]reflect.Value{reflect.ValueOf(z)})
	}
}

// benchmarkWCDMNdistance is a helper function to be called by specific benchmarkWCDMs
func benchmarkWCDMNdistance(n int, distFunc string, b *testing.B) {
	cos := WCDM{H0: 70, Om0: 0.2, Ol0: 0.7, W0: -1.2}
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

func BenchmarkWCDMComovingDistance(b *testing.B) {
	benchmarkWCDMDistance("ComovingDistance", b)
}

func BenchmarkWCDMComovingTransverseDistance(b *testing.B) {
	benchmarkWCDMDistance("ComovingTransverseDistance", b)
}

func BenchmarkWCDMLuminosityDistance(b *testing.B) {
	benchmarkWCDMDistance("LuminosityDistance", b)
}

func BenchmarkWCDMLookbackTime(b *testing.B) {
	benchmarkWCDMDistance("LookbackTime", b)
}

func BenchmarkWCDMNComovingDistance(b *testing.B) {
	benchmarkWCDMNdistance(10000, "ComovingDistance", b)
}

func BenchmarkWCDMNLuminosityDistance(b *testing.B) {
	benchmarkWCDMNdistance(10000, "LuminosityDistance", b)
}

func BenchmarkWCDMNE(b *testing.B) {
	benchmarkWCDMNdistance(10000, "E", b)
}

func BenchmarkWCDMComovingDistanceOM(b *testing.B) {
	benchmarkWCDMDistanceOM("ComovingDistance", b)
}

func BenchmarkWCDMLookbackTimeOM(b *testing.B) {
	benchmarkWCDMDistanceOM("LookbackTime", b)
}

func BenchmarkWCDMComovingDistanceFlat(b *testing.B) {
	benchmarkWCDMDistanceFlat("ComovingDistance", b)
}

func BenchmarkWCDMComovingDistanceLCDM(b *testing.B) {
	benchmarkWCDMDistanceLCDM("ComovingDistance", b)
}
