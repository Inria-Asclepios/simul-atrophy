/*=========================================================================
 *
 *  Copyright Insight Software Consortium
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *         http://www.apache.org/licenses/LICENSE-2.0.txt
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 *
 *=========================================================================*/
#ifndef __itkGradientNDAnisotropicDiffusionWithMaskFunction_h
#define __itkGradientNDAnisotropicDiffusionWithMaskFunction_h

#include "itkScalarAnisotropicDiffusionWithMaskFunction.h"
#include "itkNeighborhoodAlgorithm.h"
#include "itkNeighborhoodInnerProduct.h"
#include "itkDerivativeOperator.h"

namespace itk
{
/** \class GradientNDAnisotropicDiffusionWithMaskFunction
 *
 * This class implements an N-dimensional version of the classic Perona-Malik
 * anisotropic diffusion equation for scalar-valued images.  See
 * itkAnisotropicDiffusionWithMaskFunction for an overview of the anisotropic diffusion
 * framework and equation.
 *
 * \par
 * The conductance term for this implementation is chosen as a function of the
 * gradient magnitude of the image at each point, reducing the strength of
 * diffusion at edge pixels.
 *
 * \f[C(\mathbf{x}) = e^{-(\frac{\parallel \nabla U(\mathbf{x}) \parallel}{K})^2}\f].
 *
 * \par
 * The numerical implementation of this equation is similar to that described
 * in the Perona-Malik paper below, but uses a more robust technique
 * for gradient magnitude estimation and has been generalized to N-dimensions.
 *
 * \par References
 * Pietro Perona and Jalhandra Malik, ``Scale-space and edge detection using
 * anisotropic diffusion,'' IEEE Transactions on Pattern Analysis Machine
 * Intelligence, vol. 12, pp. 629-639, 1990.
 *
 * \sa AnisotropicDiffusionWithMaskFunction
 * \sa VectorAnisotropicDiffusionFunction
 * \sa VectorGradientAnisotropicDiffusionFunction
 * \sa CurvatureNDAnisotropicDiffusionFunction
 * \ingroup FiniteDifferenceFunctions
 * \ingroup ImageEnhancement
 * \ingroup ITKAnisotropicSmoothing
 */
template< typename TImage >
class GradientNDAnisotropicDiffusionWithMaskFunction:
  public ScalarAnisotropicDiffusionWithMaskFunction< TImage >
{
public:
  /** Standard class typedefs. */
  typedef GradientNDAnisotropicDiffusionWithMaskFunction       Self;
  typedef ScalarAnisotropicDiffusionWithMaskFunction< TImage > Superclass;
  typedef SmartPointer< Self >                         Pointer;
  typedef SmartPointer< const Self >                   ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods) */
  itkTypeMacro(GradientNDAnisotropicDiffusionWithMaskFunction,
               ScalarAnisotropicDiffusionWithMaskFunction);

  /** Inherit some parameters from the superclass type. */
  typedef typename Superclass::ImageType        ImageType;
  typedef typename Superclass::PixelType        PixelType;
  typedef typename Superclass::PixelRealType    PixelRealType;
  typedef typename Superclass::TimeStepType     TimeStepType;
  typedef typename Superclass::RadiusType       RadiusType;
  typedef typename Superclass::NeighborhoodType NeighborhoodType;
  typedef typename Superclass::FloatOffsetType  FloatOffsetType;

  typedef SizeValueType NeighborhoodSizeValueType;

  /** Inherit some parameters from the superclass type. */
  itkStaticConstMacro(ImageDimension, unsigned int, Superclass::ImageDimension);

  /** Compute the equation value. */
  virtual PixelType ComputeUpdate(const NeighborhoodType & neighborhood,
                                  void *globalData,
                                  const FloatOffsetType & offset = FloatOffsetType(0.0)
                                  );

  /** This method is called prior to each iteration of the solver. */
  virtual void InitializeIteration()
  {
    m_K = static_cast< PixelType >( this->GetAverageGradientMagnitudeSquared()
                                    * this->GetConductanceParameter() * this->GetConductanceParameter() * -2.0f );
    //BishChange:
//    m_K = -0.1;
  }

//  /** This method sets the conductance image **/
//  void SetConductanceImage(typename ImageType::Pointer conductanceImage) {
//      m_ConductanceImage = conductanceImage;
//      m_IsConductanceImageSet = true;
//  }

protected:
  GradientNDAnisotropicDiffusionWithMaskFunction();
  ~GradientNDAnisotropicDiffusionWithMaskFunction() {}

  /** Conductance Image, Conductance Image Type is same as the input image type (For now!) **/
//  typename TImage::Pointer          m_ConductanceImage;
//  bool                              m_IsConductanceImageSet;


  /** Inner product function. */
  NeighborhoodInnerProduct< ImageType > m_InnerProduct;

  /** Slices for the ND neighborhood. */
  std::slice x_slice[ImageDimension];
  std::slice xa_slice[ImageDimension][ImageDimension];
  std::slice xd_slice[ImageDimension][ImageDimension];

  /** Derivative operator. */
  DerivativeOperator< PixelType, itkGetStaticConstMacro(ImageDimension) > dx_op;

  /** Modified global average gradient magnitude term. */
  PixelType m_K;

  NeighborhoodSizeValueType m_Center;
  NeighborhoodSizeValueType m_Stride[ImageDimension];

  static double m_MIN_NORM;

private:
  GradientNDAnisotropicDiffusionWithMaskFunction(const Self &); //purposely not
                                                        // implemented
  void operator=(const Self &);                         //purposely not

  // implemented
};
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkGradientNDAnisotropicDiffusionWithMaskFunction.hxx"
#endif

#endif
