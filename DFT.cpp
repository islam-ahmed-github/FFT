

/*
 * DFT.cpp
 *
 *  Created on: Dec 15, 2018
 *      Author: Islam Ahmed
 */




#include "DFT.h"



ComplexMath::ComplexMath()
{
	//ctor
}

ComplexMath::~ComplexMath()
{
	//dtor
}



ComplexMath::Complex ComplexMath::Sum(const ComplexMath::Complex *a, const ComplexMath::Complex *b)
{
	ComplexMath::Complex temp(0, 0);
	temp.Real = (a->Real + b->Real);
	temp.Imaginary = (a->Imaginary + b->Imaginary);
	return temp;
}


ComplexMath::Complex ComplexMath::Sub(const ComplexMath::Complex *a, const ComplexMath::Complex *b)
{
	ComplexMath::Complex temp(0, 0);
	temp.Real = (a->Real - b->Real);
	temp.Imaginary = (a->Imaginary - b->Imaginary);
	return temp;
}


ComplexMath::Complex ComplexMath::Mul(const ComplexMath::Complex *a, const ComplexMath::Complex *b)
{
	ComplexMath::Complex temp(0, 0);
	temp.Real = ((a->Real * b->Real) - (a->Imaginary * b->Imaginary));
	temp.Imaginary = ((a->Real * b->Imaginary) + (a->Imaginary * b->Real));
	return temp;
}


double ComplexMath::GetMagnitude(const ComplexMath::Complex *a)
{
	double temp = 0;
	temp = (a->Real * a->Real) + (a->Imaginary * a->Imaginary);
	return  std::sqrt(temp);
}





DFT::DFT()
{

}

DFT::~DFT()
{

}


void DFT::SetUpTwiddlesLookupTable(const uint32_t FFTSize , struct Complex *OuputTwiddles)
{


	//uint64_t  f = 0, s = 0;
	//double value = 0;
	//double Cterm = ((2 * PI) / FFTSize);

	//
	//for (f = 0; f < FFTSize; f++)
	//{4
	//	for (s = 0; s < FFTSize; s++)
	//	{
	//		value = (Cterm * f * s);
	//		OuputTwiddles[MACRO_GetArrayIndex(f,s,FFTSize)]->    //= ComplexMath::Complex(cos(value), (-1 * sin(value)));
	//	}
	//}

}



void DFT::ApplyHanningWindow(const double* Src, double* Dest, const uint32_t Len)
{
	uint32_t  f = 0;
	double Cterm = ((2 * PI) / Len);
	for (f = 0; f < Len; f++) Dest[f] = (Src[f])*(0.5 - (0.5 * cos(Cterm*f)));
}



void DFT::ApplyHammingWindow(const double* Src, double* Dest, const uint32_t Len)
{
	uint64_t  f = 0;
	double Cterm = ((2 * PI) / Len);
	for (f = 0; f < Len; f++) Dest[f] = (Src[f])*(0.54 - (0.46 * cos(Cterm*f)));
}





void DFT::GetMagnitudes(const ComplexMath::Complex* Src, uint32_t Len, double* Dest)
{

	uint32_t  i = 0;
	ComplexMath c = ComplexMath();
	for (i = 0; i < Len; i++)  Dest[i] = c.GetMagnitude(&Src[i]);
}



void DFT::GetPowerSpectral(const double* SrcFrame,  double* Dest, const uint32_t Len)
{

	uint32_t  t = 0;
	for (t = 0; t <Len; t++) Dest[t] = ((SrcFrame[t] * SrcFrame[t]) / Len);

}





void DFT::Forward(const double *Input, uint32_t InputLength, double *Output, const bool HalfLength)
{
	uint64_t f = 0, s = 0, Len = 0;
	ComplexMath Cobj = ComplexMath();
	ComplexMath::Complex c1(0, 0);
	ComplexMath::Complex c2(0, 0);
	ComplexMath::Complex  P(0, 0);
	double Cterm = ((2 * PI) / InputLength);
	if (HalfLength) Len = ((InputLength / 2) + 1);
	else Len = InputLength;
	ComplexMath::Complex output(0, 0);
	for (f = 0; f < Len; f++)
	{
		output.Real = 0;
		output.Imaginary = 0;
		for (s = 0; s < InputLength; s++)
		{
			c1 = ComplexMath::Complex(Input[s], 0);
			c2 = ComplexMath::Complex(cos(Cterm*f*s), (-1 * sin(Cterm*f*s))); // twiddle
			P = Cobj.Mul(&c1, &c2);
			output = Cobj.Sum(&output, &P);
		}
		Output[f] = Cobj.GetMagnitude(&output);

	}

}//end function




void DFT::Backward(const double *Input, const uint32_t InputLength, double *Output, const bool HalfLength)
{
	uint32_t f = 0, s = 0, Len = 0;
	double Cterm = ((2 * PI) / InputLength);

	ComplexMath Cobj = ComplexMath();
	ComplexMath::Complex c1(0, 0);
	ComplexMath::Complex c2(0, 0);
	ComplexMath::Complex  P(0, 0);

	if (HalfLength) Len = ((InputLength / 2) + 1);
	else Len = InputLength;
	ComplexMath::Complex output(0, 0);
	for (f = 0; f < Len; f++)
	{
		output.Real = 0;
		output.Imaginary = 0;
		for (s = 0; s < InputLength; s++)
		{
			c1 = ComplexMath::Complex(Input[s], 0);
			c2 = ComplexMath::Complex(cos(Cterm*f*s), (sin(Cterm*f*s))); // twiddle
			P = Cobj.Mul(&c1, &c2);
			output = Cobj.Sum(&output, &P);
		}
		Output[f] = Cobj.GetMagnitude(&output);

	}

}//end function




