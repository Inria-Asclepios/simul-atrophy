#include "InverseDisplacementImageFilter.h"

//** GAUSSIANTRANSFORMNUMERICALINVERSEIMAGEFILTER **//

template<typename TDisplacementField >
InverseDisplacementImageFilter<TDisplacementField >::InverseDisplacementImageFilter(): m_MaximumNumberOfIterations(100), m_ErrorTolerance(1e-2),
    m_NumberOfErrorToleranceFailures(0)
{
}

/* Does the real work */
template< typename TDisplacementField >
void InverseDisplacementImageFilter<TDisplacementField >::ThreadedGenerateData(const OutputImageRegionType & outputRegionForThread, itk::ThreadIdType threadId)
{
	// In the fixed point iteration, we will need to access non-grid points.
	// Currently, the best interpolator that is supported by itk for vector
	// images is the linear interpolater.
	typename FieldInterpolatorType::Pointer vectorInterpolator =
			FieldInterpolatorType::New();
    vectorInterpolator->SetInputImage(this->GetInput());

	DisplacementFieldPointerType output = this->GetOutput();
	RealValueType sqr_tol = m_ErrorTolerance*m_ErrorTolerance;

	itk::ImageRegionIteratorWithIndex< DisplacementFieldType > it_dis(output, outputRegionForThread);
	for(it_dis.GoToBegin(); !it_dis.IsAtEnd(); ++it_dis){
		PointType x;
		output->TransformIndexToPhysicalPoint(it_dis.GetIndex(), x);

        VectorType v = -vectorInterpolator->Evaluate(x); // change here m_Transform->Direct(x) by the interpolated value of the input displacement field at x
		VectorType eps = v;
		unsigned int i = 2;
		
		while( (i<m_MaximumNumberOfIterations) && (eps.GetSquaredNorm()>sqr_tol) )
		{
            PointType X = x + v;
            v = -vectorInterpolator->Evaluate(X);
            PointType y = X - v; // change here m_Transform->Direct(z) by the interpolated value of the input displacement field at z

			eps = y - x;
            //v = z - y;
            ++i;
		}

		it_dis.Set( v );

        // Find how many did not converge.
		if( eps.GetSquaredNorm()>sqr_tol ){
            ++m_NumberOfErrorToleranceFailures;
		}
	}
//    std::cout << "Tolerance not reached for "<<m_NumberOfErrorToleranceFailures<<" pixels" << std::endl;
}
