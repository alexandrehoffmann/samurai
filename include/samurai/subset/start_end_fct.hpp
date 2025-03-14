// Copyright 2018-2025 the samurai's authors
// SPDX-License-Identifier:  BSD-3-Clause

#pragma once

#include "utils.hpp"

namespace samurai
{

    inline auto default_function()
    {
        return [](auto, auto i)
        {
            return i;
        };
    }

    template <std::size_t dim>
    struct start_end_function
    {
        auto& operator()(std::size_t level, std::size_t min_level, std::size_t max_level)
        {
            m_level = level;
            m_shift = static_cast<int>(max_level) - static_cast<int>(min_level);
            return *this;
        }

        template <std::size_t, class Func>
        inline auto start(const Func& f) const
        {
            auto new_f = [&, f](auto, auto i)
            {
                return f(m_level, (i >> m_shift) << m_shift);
            };
            return new_f;
        }

        template <std::size_t, class Func>
        inline auto end(const Func& f) const
        {
            auto new_f = [&, f](auto, auto i)
            {
                return f(m_level, (((i - 1) >> m_shift) + 1) << m_shift);
            };
            return new_f;
        }

        template <std::size_t, class Func>
        inline auto goback(const Func& f) const
        {
            auto new_f = [&, f](auto level, auto i)
            {
                i = start_shift(f(m_level, i), static_cast<int>(level) - static_cast<int>(m_level));
                return i;
            };
            return new_f;
        }

        std::size_t m_level;
        int m_shift;
    };

    template <std::size_t dim>
    struct start_end_translate_function
    {
        using container_t = xt::xtensor_fixed<int, xt::xshape<dim>>;

        explicit start_end_translate_function(const container_t& t)
            : m_level(0)
            , m_min_level(0)
            , m_max_level(0)
            , m_t(t)
        {
        }

        auto& operator()(auto level, auto min_level, auto max_level)
        {
            m_level     = level;
            m_min_level = min_level;
            m_max_level = max_level;
            return *this;
        }

        template <std::size_t d, class Func>
        inline auto start(const Func& f) const
        {
            auto new_f = [&, f](auto level, auto i)
            {
                int max2curr = static_cast<int>(m_max_level) - static_cast<int>(level);
                int curr2min = static_cast<int>(level) - static_cast<int>(m_min_level);
                int min2max  = static_cast<int>(m_max_level) - static_cast<int>(m_min_level);

                return f(m_level, (((i >> max2curr) + m_t[d - 1]) >> curr2min) << min2max);
            };
            return new_f;
        }

        template <std::size_t d, class Func>
        inline auto end(const Func& f) const
        {
            auto new_f = [&, f](auto level, auto i)
            {
                int max2curr = static_cast<int>(m_max_level) - static_cast<int>(level);
                int curr2min = static_cast<int>(level) - static_cast<int>(m_min_level);
                int min2max  = static_cast<int>(m_max_level) - static_cast<int>(m_min_level);

                return f(m_level, (((((i - 1) >> max2curr) + m_t[d - 1]) >> curr2min) + 1) << min2max);
            };
            return new_f;
        }

        template <std::size_t d, class Func>
        inline auto goback(const Func& f) const
        {
            auto new_f = [&, f](auto level, auto i)
            {
                i = start_shift(f(m_level, i), static_cast<int>(level) - static_cast<int>(m_level)) - m_t[d - 1];
                return i;
            };
            return new_f;
        }

        std::size_t m_level;
        std::size_t m_min_level;
        std::size_t m_max_level;
        xt::xtensor_fixed<int, xt::xshape<dim>> m_t;
    };

    template <std::size_t dim>
    struct start_end_contraction_function
    {
        explicit start_end_contraction_function(int c)
            : m_level(0)
            , m_min_level(0)
            , m_max_level(0)
            , m_c(c)
        {
        }

        auto& operator()(auto level, auto min_level, auto max_level)
        {
            m_level     = level;
            m_min_level = min_level;
            m_max_level = max_level;
            return *this;
        }

        template <std::size_t d, class Func>
        inline auto start(const Func& f) const
        {
            auto new_f = [&, f](auto level, auto i)
            {
                int max2curr = static_cast<int>(m_max_level) - static_cast<int>(level);
                int curr2min = static_cast<int>(level) - static_cast<int>(m_min_level);
                int min2max  = static_cast<int>(m_max_level) - static_cast<int>(m_min_level);

                return f(m_level, (((i >> max2curr) - m_c) >> curr2min) << min2max);
            };
            return new_f;
        }

        template <std::size_t d, class Func>
        inline auto end(const Func& f) const
        {
            auto new_f = [&, f](auto level, auto i)
            {
                int max2curr = static_cast<int>(m_max_level) - static_cast<int>(level);
                int curr2min = static_cast<int>(level) - static_cast<int>(m_min_level);
                int min2max  = static_cast<int>(m_max_level) - static_cast<int>(m_min_level);

                return f(m_level, (((((i - 1) >> max2curr) - m_c) >> curr2min) + 1) << min2max);
            };
            return new_f;
        }

        template <std::size_t d, class Func>
        inline auto goback(const Func& f) const
        {
            auto new_f = [&, f](auto level, auto i)
            {
                i = start_shift(f(m_level, i), static_cast<int>(level) - static_cast<int>(m_level)) + m_c;
                return i;
            };
            return new_f;
        }

        std::size_t m_level;
        std::size_t m_min_level;
        std::size_t m_max_level;
        int m_c;
    };
}
