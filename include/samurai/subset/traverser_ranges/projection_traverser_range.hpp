// Copyright 2018-2025 the samurai's authors
// SPDX-License-Identifier:  BSD-3-Clause

#pragma once

#include "set_traverser_range_base.hpp"

namespace samurai
{

    template <class SetTraverserRange>
    class ProjectionTraverserRange;

    template <class SetTraverserRange>
    struct SetTraverserRangeTraits<ProjectionTraverserRange<SetTraverserRange>>
    {
    };

} // namespace samurai
