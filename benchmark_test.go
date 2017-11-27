package cosmo

import (
	"testing"
)

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
