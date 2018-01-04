package cosmo

import (
	"reflect"
	"testing"
)

func benchmarkFlatLCDMEN(n int, b *testing.B) {
	cos := FlatLCDM{Om0: 0.27, H0: 70}

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

func BenchmarkFlatLCDMEN(b *testing.B) {
	benchmarkFlatLCDMEN(10000, b)
}

func BenchmarkFlatLCDMENdistance(b *testing.B) {
	benchmarkFlatLCDMNdistance(10000, "E", b)
}

func BenchmarkFlatLCDME(b *testing.B) {
	cos := FlatLCDM{Om0: 0.27, H0: 70}
	z := 1.0
	for i := 0; i < b.N; i++ {
		cos.E(z)
	}
}

func BenchmarkFlatLCDMEinv(b *testing.B) {
	cos := FlatLCDM{Om0: 0.27, H0: 70}
	z := 1.0
	for i := 0; i < b.N; i++ {
		cos.Einv(z)
	}
}

// benchmarkFlatLCDMDistance is a helper function to be called by specific benchmarkFlatLCDMs
func benchmarkFlatLCDMDistance(distFunc string, b *testing.B) {
	cos := FlatLCDM{Om0: 0.27, H0: 70}
	z := 1.0

	funcToTest := reflect.ValueOf(&cos).MethodByName(distFunc)
	for i := 0; i < b.N; i++ {
		funcToTest.Call([]reflect.Value{reflect.ValueOf(z)})
	}
}

// benchmarkFlatLCDMNdistance is a helper function to be called by specific benchmarkFlatLCDMs
func benchmarkFlatLCDMNdistance(n int, distFunc string, b *testing.B) {
	cos := FlatLCDM{Om0: 0.27, H0: 70}
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

func BenchmarkFlatLCDMComovingDistance(b *testing.B) {
	benchmarkFlatLCDMDistance("ComovingDistance", b)
}

func BenchmarkFlatLCDMComovingTransverseDistance(b *testing.B) {
	benchmarkFlatLCDMDistance("ComovingTransverseDistance", b)
}

func BenchmarkFlatLCDMLuminosityDistance(b *testing.B) {
	benchmarkFlatLCDMDistance("LuminosityDistance", b)
}

func BenchmarkFlatLCDMLookbackTime(b *testing.B) {
	benchmarkFlatLCDMDistance("LookbackTime", b)
}

func BenchmarkFlatLCDMNComovingDistance(b *testing.B) {
	benchmarkFlatLCDMNdistance(10000, "ComovingDistance", b)
}

func BenchmarkFlatLCDMNLuminosityDistance(b *testing.B) {
	benchmarkFlatLCDMNdistance(10000, "LuminosityDistance", b)
}

func BenchmarkFlatLCDMNE(b *testing.B) {
	benchmarkFlatLCDMNdistance(10000, "E", b)
}
