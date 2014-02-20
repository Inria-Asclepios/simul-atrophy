#ifndef INVERSE_DISPLACEMENT_FILTER_H
#define INVERSE_DISPLACEMENT_FILTER_H

/*
    class InverseDisplacementImageFilter
        Numerical "inverse" for (possibly non-invertible) displacement fields. Uses a simple fixed point scheme that should do well for relatively smooth fields:
			v_0 = 0
			v_{i+1} = -u( x + v_i )
		The input image is simply here to define the domain and resolution of interest over which to sample the inverse field.

		itk::Vector<T,d> pixel type expected for the displacement field.
*/
#include "itkImageToImageFilter.h"
//#include "itkVectorLinearInterpolateImageFunction.h"
#include <itkVectorLinearInterpolateNearestNeighborExtrapolateImageFunction.h>
#include "itkImageRegionIterator.h"
#include "itkImageRegionIteratorWithIndex.h"

template<typename TDisplacementField>
class InverseDisplacementImageFilter: public itk::ImageToImageFilter< TDisplacementField, TDisplacementField >
{
	public:
	/* Standard class typedefs. */
    typedef InverseDisplacementImageFilter										Self;
	typedef itk::ImageToImageFilter< TDisplacementField, TDisplacementField >		Superclass;
	typedef itk::SmartPointer< Self >																		Pointer;
	typedef itk::SmartPointer< const Self >																ConstPointer;

	/* ImageDimension constants */
	itkStaticConstMacro(Dimension, unsigned int, TDisplacementField::ImageDimension);

	/* Image typedefs */
	typedef TDisplacementField											DisplacementFieldType;
	typedef typename DisplacementFieldType::Pointer			DisplacementFieldPointerType;

	typedef double											SpacingValueType;
	typedef double											RealValueType;

	/* Container typedefs */
	typedef typename DisplacementFieldType::PointType		PointType;
	typedef typename DisplacementFieldType::PixelType		VectorType;

	/* Interpolator typedef */
    typedef itk::VectorLinearInterpolateNearestNeighborExtrapolateImageFunction<DisplacementFieldType,double> FieldInterpolatorType;
    typedef typename FieldInterpolatorType::Pointer                  FieldInterpolatorPointer;
  	typedef typename FieldInterpolatorType::OutputType               FieldInterpolatorOutputType;

    typedef typename Superclass::OutputImageRegionType OutputImageRegionType;

	/* Method for creation through the object factory. */
	itkNewMacro(Self);
 
	/* Run-time type information (and related methods). */
    itkTypeMacro(InverseDisplacementImageFilter, itk::ImageToImageFilter);

	/* Set\get methods */
	itkSetMacro(MaximumNumberOfIterations, unsigned int);
	itkSetMacro(ErrorTolerance, RealValueType);

	itkGetMacro(MaximumNumberOfIterations, unsigned int);
	itkGetMacro(ErrorTolerance, RealValueType);

    itkGetMacro(NumberOfErrorToleranceFailures, unsigned int);

protected:
    InverseDisplacementImageFilter();
    ~InverseDisplacementImageFilter() {}

	/** Does the real work. */
    virtual void ThreadedGenerateData(const OutputImageRegionType & outputRegionForThread, itk::ThreadIdType threadId);

private:
    InverseDisplacementImageFilter(const Self &);	//purposely not implemented
	void operator=(const Self &);								//purposely not implemented

	unsigned int m_MaximumNumberOfIterations;					// hard cut regardless of convergence in the fixed point scheme
	RealValueType m_ErrorTolerance;								// if we go below this tolerance, convergence is reached.
    unsigned int m_NumberOfErrorToleranceFailures;              // total number of voxels for which convergence was not reached.
};

#include "InverseDisplacementImageFilter.txx"

#endif
