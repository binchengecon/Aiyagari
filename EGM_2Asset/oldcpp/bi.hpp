/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2003 Ferdinando Ametrano
 Copyright (C) 2005 Gary Kennedy

 This file is part of QuantLib, a free-software/open-source library
 for financial quantitative analysts and developers - http://quantlib.org/

 QuantLib is free software: you can redistribute it and/or modify it
 under the terms of the QuantLib license.  You should have received a
 copy of the license along with this program; if not, please email
 <quantlib-dev@lists.sf.net>. The license is also available online at
 <http://quantlib.org/license.shtml>.

 This program is distributed in the hope that it will be useful, but WITHOUT
 ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 FOR A PARTICULAR PURPOSE.  See the license for more details.
*/

#ifndef quantlib_bivariatenormal_distribution_hpp
#define quantlib_bivariatenormal_distribution_hpp

#include <ql/math/distributions/normaldistribution.hpp>

namespace QuantLib
{

    class BivariateCumulativeNormalDistributionDr78
    {
    public:
        BivariateCumulativeNormalDistributionDr78(Real rho);
        // function
        Real operator()(Real a, Real b) const;

    private:
        Real rho_, rho2_;
        static const Real x_[], y_[];
    };

    class BivariateCumulativeNormalDistributionWe04DP
    {
    public:
        BivariateCumulativeNormalDistributionWe04DP(Real rho);
        // function
        Real operator()(Real a, Real b) const;

    private:
        Real correlation_;
        CumulativeNormalDistribution cumnorm_;
    };

    typedef BivariateCumulativeNormalDistributionWe04DP
        BivariateCumulativeNormalDistribution;

}

#endif