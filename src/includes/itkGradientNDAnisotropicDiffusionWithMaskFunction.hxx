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
#ifndef __itkGradientNDAnisotropicDiffusionWithMaskFunction_hxx
#define __itkGradientNDAnisotropicDiffusionWithMaskFunction_hxx

#include "itkNumericTraits.h"
#include "itkGradientNDAnisotropicDiffusionWithMaskFunction.h"

namespace itk
{
template< typename TImage >
double GradientNDAnisotropicDiffusionWithMaskFunction< TImage >
::m_MIN_NORM = 1.0e-10;

template< typename TImage >
GradientNDAnisotropicDiffusionWithMaskFunction< TImage >
::GradientNDAnisotropicDiffusionWithMaskFunction() :
    m_K(0.0)
{
    //    m_IsConductanceImageSet = false;
    unsigned int i, j;
    RadiusType   r;

    for ( i = 0; i < ImageDimension; ++i )
    {
        r[i] = 1;
    }
    this->SetRadius(r);

    // Dummy neighborhood used to set up the slices.
    Neighborhood< PixelType, ImageDimension > it;
    it.SetRadius(r);

    // Slice the neighborhood
    m_Center =  it.Size() / 2;

    for ( i = 0; i < ImageDimension; ++i )
    {
        m_Stride[i] = it.GetStride(i);
    }

    for ( i = 0; i < ImageDimension; ++i )
    {
        x_slice[i]  = std::slice(m_Center - m_Stride[i], 3, m_Stride[i]);
    }

    for ( i = 0; i < ImageDimension; ++i )
    {
        for ( j = 0; j < ImageDimension; ++j )
        {
            // For taking derivatives in the i direction that are offset one
            // pixel in the j direction.
            xa_slice[i][j] =
                    std::slice( ( m_Center + m_Stride[j] ) - m_Stride[i], 3, m_Stride[i] );
            xd_slice[i][j] =
                    std::slice( ( m_Center - m_Stride[j] ) - m_Stride[i], 3, m_Stride[i] );
        }
    }

    // Allocate the derivative operator.
    dx_op.SetDirection(0);  // Not relevant, will be applied in a slice-based
    // fashion.
    dx_op.SetOrder(1);
    dx_op.CreateDirectional();
}

template< typename TImage >
typename GradientNDAnisotropicDiffusionWithMaskFunction< TImage >::PixelType
GradientNDAnisotropicDiffusionWithMaskFunction< TImage >
::ComputeUpdate(const NeighborhoodType & it, void *,
                const FloatOffsetType &)
{
    unsigned int i, j;

    double accum;
    double accum_d;
    double Cx;
    double Cxd;

    // PixelType is scalar in this context
    PixelRealType delta;
    PixelRealType dx_forward;
    PixelRealType dx_backward;
    //    PixelRealType dx[ImageDimension];
    PixelRealType dx_aug;
    PixelRealType dx_dim;

    // Iterator for the conductance image

    NeighborhoodType          it_conductance(it.GetRadius(), this->m_ConductanceImage,it.GetRegion());
    it_conductance.SetLocation(it.GetIndex());

    //      it_conductance.SetBoundaryCondition(it.GetBoundaryCondition());
    //    // Derivative terms for conductance Image
    PixelRealType dx_forward_c;
    PixelRealType dx_backward_c;
    PixelRealType dx_c[ImageDimension];
    //    PixelRealType dx_aug_c;
    //    PixelRealType dx_dim_c;

    delta = NumericTraits< PixelRealType >::ZeroValue();

    // Calculate the centralized derivatives for each dimension.
    for ( i = 0; i < ImageDimension; i++ )
    {
        //        dx[i]  =  ( it.GetPixel(m_Center + m_Stride[i]) - it.GetPixel(m_Center - m_Stride[i]) ) / 2.0f;
        //        dx[i] *= this->m_ScaleCoefficients[i];
        dx_c[i]  =  ( it_conductance.GetPixel(m_Center + m_Stride[i]) - it_conductance.GetPixel(m_Center - m_Stride[i]) ) / 2.0f;
        dx_c[i] *= this->m_ScaleCoefficients[i];
        //        if(fabs(dx[i]) > 1e-6) std::cout<<"dx["<<i<<"]="<<dx[i]<<"\t";    OK
        //        if(fabs(dx_c[i]) > 1e-6) std::cout<<"dx_c["<<i<<"]="<<dx_c[i]<<"\n";  OK

    }

    for ( i = 0; i < ImageDimension; i++ )
    {
        // ``Half'' directional derivatives
        dx_forward = it.GetPixel(m_Center + m_Stride[i])
                - it.GetPixel(m_Center);
        dx_forward *= this->m_ScaleCoefficients[i];
        dx_backward =  it.GetPixel(m_Center)
                - it.GetPixel(m_Center - m_Stride[i]);
        dx_backward *= this->m_ScaleCoefficients[i];

        dx_forward_c = it_conductance.GetPixel(m_Center + m_Stride[i])
                - it_conductance.GetPixel(m_Center);
        dx_forward_c *= this->m_ScaleCoefficients[i];
        dx_backward_c =  it_conductance.GetPixel(m_Center)
                - it_conductance.GetPixel(m_Center - m_Stride[i]);
        dx_backward_c *= this->m_ScaleCoefficients[i];

        // Calculate the conductance terms.  Conductance varies with each
        // dimension because the gradient magnitude approximation is different
        // along each  dimension.
        accum   = 0.0;
        accum_d = 0.0;
        for ( j = 0; j < ImageDimension; j++ )
        {
            if ( j != i )
            {
                dx_aug = ( it_conductance.GetPixel(m_Center + m_Stride[i] + m_Stride[j])
                           - it_conductance.GetPixel(m_Center + m_Stride[i] - m_Stride[j]) ) / 2.0f;
                dx_aug *= this->m_ScaleCoefficients[j];
                dx_dim = ( it_conductance.GetPixel(m_Center - m_Stride[i] + m_Stride[j])
                           - it_conductance.GetPixel(m_Center - m_Stride[i] - m_Stride[j]) ) / 2.0f;
                dx_dim *= this->m_ScaleCoefficients[j];
                accum += 0.25f * vnl_math_sqr(dx_c[j] + dx_aug);
                accum_d += 0.25f * vnl_math_sqr(dx_c[j] + dx_dim);
            }
        }
        if ( m_K == 0.0 )
        {
            Cx = 0.0;
            Cxd = 0.0;
        }
        else
        {
            Cx = std::exp( ( vnl_math_sqr(dx_forward_c) + accum )  / m_K );
            Cxd = std::exp( ( vnl_math_sqr(dx_backward_c) + accum_d ) / m_K );
        }
        //Calcualtion of conductance done, now use them to modify the derivates
        // Conductance modified first order derivatives.
        dx_forward  = dx_forward * Cx;
        dx_backward = dx_backward * Cxd;

        // Conductance modified second order derivative.
        delta += dx_forward - dx_backward;
    }

    return static_cast< PixelType >( delta );
}
} // end namespace itk

#endif
