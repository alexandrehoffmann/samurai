// Copyright 2018-2025 the samurai's authors
// SPDX-License-Identifier:  BSD-3-Clause

#pragma once

#include "set_traverser_range_base.hpp"

namespace samurai
{

    template <class SetTraverserRange>
    class TranslationTraverserRange;

    template <class SetTraverserRange>
    struct SetTraverserRangeTraits<TranslationTraverserRange<SetTraverserRange>>
    {
        static_assert(IsSetTraverserRange<SetTraverserRange>::value);

        using ChildIterator = typename SetTraverserRange::Iterator;
        using Translation   = typename ChildIterator::value_type::value_type;

        class Iterator
        {
          public:

            using iterator_category = std::forward_iterator_tag;
            using difference_type   = std::ptrdiff_t;
            using value_type        = TranslationTraverser<typename SetTraverserRange::Element>;
            using reference         = value_type;

            Iterator(const ChildIterator childIterator, const Translation& translation)
                : m_childIterator(childIterator)
                , m_translation(translation)
            {
            }

            reference operator*() const
            {
                return reference(*m_childIterator, m_translation);
            }

            Iterator& operator++()
            {
                ++m_childIterator;
                return *this;
            }

            Iterator operator++(int)
            {
                Iterator tmp = *this;
                ++(*this);
                return tmp;
            }

            friend bool operator==(const Iterator& a, const Iterator& b)
            {
                return a.m_childIterator == b.m_childIterator;
            };

            friend bool operator!=(const Iterator& a, const Iterator& b)
            {
                return a.m_childIterator != b.m_childIterator;
            };

          private:

            ChildIterator m_childIterator;
            Translation m_translation;
        };
    };

    template <class SetTraverserRange>
    class TranslationTraverserRange : public SetTraverserRangeBase<TranslationTraverserRange<SetTraverserRange>>
    {
        using Self = TranslationTraverserRange<SetTraverserRange>;

      public:

        SAMURAI_SET_TRAVERSER_RANGE_TYPEDEFS

        using Translation = typename SetTraverserRangeTraits<Self>::Translation;

        TranslationTraverserRange(const SetTraverserRange& set_traverser_range, const Translation& translation)
            : m_set_traverser_range(set_traverser_range)
            , m_translation(translation)
        {
        }

        Iterator begin_impl()
        {
            return Iterator(m_set_traverser_range.begin(), m_translation);
        }

        Iterator end_impl()
        {
            return Iterator(m_set_traverser_range.end(), m_translation);
        }

      private:

        SetTraverserRange m_set_traverser_range;
        Translation m_translation;
    };

} // namespace samurai
