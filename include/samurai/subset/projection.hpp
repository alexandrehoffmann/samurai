// Copyright 2018-2025 the samurai's authors
// SPDX-License-Identifier:  BSD-3-Clause

#pragma once

#include "../samurai_config.hpp"
#include "../static_algorithm.hpp"
#include "set_base.hpp"
#include "traversers/projection_traverser.hpp"

namespace samurai
{

    template <class Set>
    class Projection;

    template <class Set>
    struct SetTraits<Projection<Set>>
    {
        static_assert(IsSet<Set>::value);

        template <std::size_t d>
        using traverser_t = ProjectionTraverser<typename Set::template traverser_range_t<d>>;

        static constexpr std::size_t dim()
        {
            return Set::dim;
        }
    };

    template <class Set>
    class Projection : public SetBase<Projection<Set>>
    {
        using Self = Projection<Set>;

      public:

        SAMURAI_SET_TYPEDEFS

        Projection(const Set& set, const std::size_t level)
            : m_set(set)
            , m_level(level)
        {
            if (m_level < m_set.level())
            {
                m_projectionType = ProjectionType::COARSEN;
                m_shift          = m_set.level() - m_level;
            }
            else
            {
                m_projectionType = ProjectionType::REFINE;
                m_shift          = m_level - m_set.level();
            }
        }

        inline std::size_t level_impl() const
        {
            return m_level;
        }

        inline bool exist_impl() const
        {
            return m_set.exist();
        }

        inline bool empty_impl() const
        {
            return m_set.empty();
        }

        template <std::size_t d>
        inline void init_get_traverser_work_impl(const std::size_t n_traversers, std::integral_constant<std::size_t, d> d_ic) const
        {
            m_set.init_get_traverser_work(n_traversers, d_ic);
        }

        template <std::size_t d>
        inline void clear_get_traverser_work_impl(std::integral_constant<std::size_t, d> d_ic) const
        {
            m_set.clear_get_traverser_work(d_ic);
        }

        template <class index_t, std::size_t d>
        inline traverser_t<d> get_traverser_impl(const index_t& index, std::integral_constant<std::size_t, d> d_ic) const
        {
            if (m_projectionType == ProjectionType::COARSEN)
            {
                if constexpr (d != Base::dim - 1)
                {
                    xt::xtensor_fixed<value_t, xt::xshape<Base::dim - 1>> index_min(index << m_shift);
                    xt::xtensor_fixed<value_t, xt::xshape<Base::dim - 1>> index_max(index << m_shift);

                    index_max[d] = ((index[d] + 1) << m_shift) - 1;

                    return traverser_t<d>(m_set.get_traverser_range(index_min, index_max, d_ic), m_projectionType, m_shift);
                }
                else
                {
                    return traverser_t<d>(m_set.get_traverser_range(index, index, d_ic), m_projectionType, m_shift);
                }
            }
            else
            {
                return traverser_t<d>(m_set.get_traverser_range(index, index, d_ic), m_projectionType, m_shift);
            }
        }

      private:

        Set m_set;
        std::size_t m_level;
        ProjectionType m_projectionType;
        std::size_t m_shift;
    };

} // namespace samurai
