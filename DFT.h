
/*
 * DFT.h
 *
 *  Created on: Dec 15, 2018
 *      Author: Islam Ahmed
 */



#pragma once

#include<stdint.h>
#include<cmath>

#define PI ((double)3.14159265358979323846)

#define MACRO_GetArrayIndex(RowIndex,ColIndex,RowLength)  (((RowIndex)*(RowLength))+(ColIndex))

class ComplexMath
{

public:
	ComplexMath();
	virtual ~ComplexMath();

	struct Complex
	{
		double Real = 0;
		double Imaginary = 0;
		Complex() {};
		Complex(double Real, double Imaginary) : Real(Real), Imaginary(Imaginary) {}

	};

	ComplexMath::Complex Sum(const ComplexMath::Complex*, const ComplexMath::Complex*);
	ComplexMath::Complex Sub(const ComplexMath::Complex*, const ComplexMath::Complex*);
	ComplexMath::Complex Mul(const ComplexMath::Complex*, const ComplexMath::Complex*);
	double GetMagnitude(const ComplexMath::Complex*);


protected:

private:

};






class DFT
{
public:
	DFT();
	virtual ~DFT();

	void SetUpTwiddlesLookupTable( uint32_t FFTSize , struct Complex * Output);

	void Forward(const double *Input, const uint32_t InputLength, double *Output, const bool HalfLength);
	void Backward(const double *Input, const uint32_t InputLength, double *Output, const bool HalfLength);
	void GetPowerSpectral(const double* SrcFrame, double* Dest, uint32_t Len);



protected:


private:

	void GetMagnitudes(const ComplexMath::Complex* Src, uint32_t Len, double* Dest);
	void ApplyHanningWindow(const double* Src, double * Dest, const uint32_t Len);
	void ApplyHammingWindow(const double* Src, double * Dest, const uint32_t Len);


};


