package cosmo

import (
	"fmt"
	"gonum.org/v1/gonum/floats"
	"runtime"
	"testing"
)

// runTests run cos_func on an array of input inputs and compares to expected expected
func runTests(testFunc func(float64) float64, inputs, expected []float64, tol float64, t *testing.T) {
	// We ask runTest to look one additional stack level down
	// to get the original caller of runTest
	stackLevel := 1
	for i, z := range inputs {
		runTest(testFunc, z, expected[i], tol, t, stackLevel)
	}
}

// runTest runs 'testFunc' on scalar 'input' and compares to 'exp'
func runTest(testFunc func(float64) float64, input, exp, tol float64, t *testing.T, stackLevel int) {
	var test_description, test_line string

	pc, file, no, ok := runtime.Caller(stackLevel + 1)
	if ok {
		details := runtime.FuncForPC(pc)
		test_description = details.Name()
		test_line = fmt.Sprintf("%s#%d", file, no)
	} else {
		test_description = "CAN'T DETERMINE TEST NAME"
		test_line = "CAN'T DETERMINE TEST LINE"
	}

	obs := testFunc(input)
	if !floats.EqualWithinAbs(obs, exp, tol) {
		t.Errorf("Failed %s at\n %s\n"+"  Expected %f, return %f",
			test_description, test_line, exp, obs)
	}
}
