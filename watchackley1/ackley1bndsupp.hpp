/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   ackley1bndsupp.hpp
 * Author: Alexander Usov
 *
 * Created on December 10, 2016, 8:16 AM
 */

#ifndef ACKLEY1BNDSUPP_HPP
#define ACKLEY1BNDSUPP_HPP

#include <math.h>
#include <box/boxutils.hpp>
#include <common/interval.hpp>
#include <cutfact/lbcutfact/boundsupp.hpp>


namespace TESTNUC {

    /**
     * Interval analisis lower bound supplier for Ackley1 function 
     */
    template <class FT> class Ackley1BoundSupplier : public NUC::BoundSupplier <FT>{
    public:

        /**
         * Constructor
         * @param n problem dimension
         */
        Ackley1BoundSupplier(int n) : mN(n) {

        }
        
        /**
         * Retrieve 
         * @param box to find the bound
         * @return bound
         */
        FT getBound(const snowgoose::Box<FT>& box) {
            FT lb = 0;
            FT ub = 0;
            for(int i=0; i < mN; i++)
            {
                FT min, max;
                snowgoose::Interval<FT>::sqr(box.mA[i], box.mB[i], &min, &max);
                snowgoose::Interval<FT>::sum(lb, ub, min, max, &lb, &ub);
            }            
            snowgoose::Interval<FT>::mult((FT)1/mN, lb, ub, &lb, &ub);
            snowgoose::Interval<FT>::sqrt(lb, ub, &lb, &ub);
            snowgoose::Interval<FT>::mult(-0.2, lb, ub, &lb, &ub);
            snowgoose::Interval<FT>::exp(lb, ub, &lb, &ub);
            snowgoose::Interval<FT>::mult(-20, lb, ub, &lb, &ub);
            FT templb = lb;
            FT tempub = ub;
            
            lb = 0;
            ub = 0;
            for(int i=0; i < mN; i++)
            {
                FT min, max;
                snowgoose::Interval<FT>::mult(2*M_PI, box.mA[i], box.mB[i], &min, &max);                
                snowgoose::Interval<FT>::cos_(min, max, &min, &max);
                snowgoose::Interval<FT>::sum(lb, ub, min, max, &lb, &ub);
            }
            snowgoose::Interval<FT>::mult((FT)1/mN, lb, ub, &lb, &ub);
            snowgoose::Interval<FT>::exp(lb, ub, &lb, &ub);
            
            snowgoose::Interval<FT>::diff(templb, tempub, lb, ub, &lb, &ub);
            snowgoose::Interval<FT>::sum(20, lb, ub, &lb, &ub);
            snowgoose::Interval<FT>::sum(M_El, lb, ub, &lb, &ub);
            return lb;
        }

    private:
        int mN;
    };

}

#endif /* ACKLEY1BNDSUPP_HPP */

