package main

import "fmt"

//6th ed. Numerical Mathematics and Computing problem 11.2.2
//Using psuedo code on pg. 470

//RK4 Solver
func RK4(t float64, x []float64, h float64, ic int, n int){
	K1 := make([]float64, ic)
	K2 := make([]float64, ic)
	K3 := make([]float64, ic)
	K4 := make([]float64, ic)
	R := make([]float64, ic)

	for i := 0; i < n; i++{
		AiryDE(t, x, K1)

		for j := 0; j < ic; j++{
			R[j] = x[j] + h/2*K1[j]
		}
		AiryDE(t+h/2, R, K2)

		for j := 0; j < ic; j++{
			R[j] = x[j] + h/2*K2[j]
		}
		AiryDE(t+h/2, R, K3)

		for j := 0; j < ic; j++{
			R[j] = x[j] + h*K3[j]
		}
		AiryDE(t+h, R, K4)

		for j := 0; j < ic; j++{
			x[j] = x[j] + h/6*(K1[j]+2*(K2[j]+K3[j])+K4[j])
		}
		t = t + h
	}
}


//Airy DE Function
func AiryDE(t float64, x []float64, k []float64){
	k[0] = x[1]
	k[1] = t*x[0]
}

func main() {
	var t float64
	var h float64
	x := make([]float64, 2)
	t = 0.0
	h = 0.045  //h = (4.5-0)/100
	x[0] = 0.355028053887817 // intial x(0)
	x[1] = -0.258819403792807 // intial x'(0)
	RK4(t, x, h, 2, 100)
	expected := 0.0003302503
	fmt.Println(x[0], " ", x[1], " ", (expected-x[0])/2)
}

//output x[0] = 0.00032966584082762685, x[1] = -0.0007190741695770324
//expected x[0] = 0.0003302503
//absolute error (0.0003302503-0.00032966584082762685)/2 = 2.92229586186586e-07