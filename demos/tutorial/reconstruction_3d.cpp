// Copyright 2018-2025 the samurai's authors
// SPDX-License-Identifier:  BSD-3-Clause

#include <chrono>
#include <filesystem>

#include <samurai/algorithm.hpp>
#include <samurai/field.hpp>
#include <samurai/io/hdf5.hpp>
#include <samurai/mr/adapt.hpp>
#include <samurai/mr/mesh.hpp>
#include <samurai/reconstruction.hpp>
#include <samurai/samurai.hpp>
#include <samurai/uniform_mesh.hpp>

namespace fs = std::filesystem;

enum class Case : int
{
    abs,
    exp,
    tanh
};

using namespace samurai::math;

template <class Mesh>
auto init(Mesh& mesh, Case& c)
{
    using mesh_id_t = typename Mesh::mesh_id_t;

    auto u = samurai::make_scalar_field<double>("u", mesh);

    samurai::for_each_interval(mesh[mesh_id_t::cells],
                               [&](std::size_t level, const auto& i, const auto& index)
                               {
                                   auto j          = index[0];
                                   auto k          = index[1];
                                   const double dx = mesh.cell_length(level);
                                   auto x          = mesh.origin_point()[0] + dx * arange<double>(i.start, i.end) + 0.5 * dx;
                                   auto y          = mesh.origin_point()[1] + j * dx + 0.5 * dx;
                                   auto z          = mesh.origin_point()[2] + k * dx + 0.5 * dx;

                                   switch (c)
                                   {
                                       case Case::abs:
                                           u(level, i, j, k) = abs(x) + std::abs(y) + std::abs(z);
                                           break;
                                       case Case::exp:
                                           u(level, i, j, k) = exp(-100 * (x * x + y * y + z * z));
                                           break;
                                       case Case::tanh:
                                           u(level, i, j, k) = tanh(50 * (abs(x) + std::abs(y) + std::abs(z))) - 1;
                                           break;
                                   }
                               });

    switch (c)
    {
        case Case::abs:
            samurai::make_bc<samurai::Dirichlet<1>>(u,
                                                    [](auto, auto, const auto& coords)
                                                    {
                                                        return std::abs(coords[0]) + std::abs(coords[1]) + std::abs(coords[2]);
                                                    });
            break;
        case Case::exp:
            samurai::make_bc<samurai::Dirichlet<1>>(
                u,
                [](auto, auto, const auto& coords)
                {
                    return std::exp(-100 * (coords[0] * coords[0] + coords[1] * coords[1] + coords[2] * coords[2]));
                });
            break;
        case Case::tanh:
            samurai::make_bc<samurai::Dirichlet<1>>(
                u,
                [](auto, auto, const auto& coords)
                {
                    return std::tanh(50 * (std::abs(coords[0]) + std::abs(coords[1]) + std::abs(coords[2]))) - 1;
                });
            break;
    }

    return u;
}

int main(int argc, char* argv[])
{
    auto& app = samurai::initialize("3d reconstruction of an adapted solution using multiresolution", argc, argv);

    constexpr size_t dim                     = 3;
    constexpr std::size_t max_stencil_width_ = 2;
    // constexpr std::size_t graduation_width_ =
    // samurai::default_config::graduation_width;
    constexpr std::size_t graduation_width_     = 2;
    constexpr std::size_t max_refinement_level_ = samurai::default_config::max_level;
    constexpr std::size_t prediction_order_     = 1;
    using MRConfig = samurai::MRConfig<dim, max_stencil_width_, graduation_width_, prediction_order_, max_refinement_level_>;

    Case test_case{Case::abs};
    const std::map<std::string, Case> map{
        {"abs",  Case::abs },
        {"exp",  Case::exp },
        {"tanh", Case::tanh}
    };

    // Adaptation parameters
    std::size_t min_level = 2;
    std::size_t max_level = 5;

    // Output parameters
    fs::path path        = fs::current_path();
    std::string filename = "reconstruction_3d";

    app.add_option("--case", test_case, "Test case")->capture_default_str()->transform(CLI::CheckedTransformer(map, CLI::ignore_case));
    app.add_option("--min-level", min_level, "Minimum level of the multiresolution")->capture_default_str()->group("Multiresolution");
    app.add_option("--max-level", max_level, "Maximum level of the multiresolution")->capture_default_str()->group("Multiresolution");
    app.add_option("--path", path, "Output path")->capture_default_str()->group("Output");
    app.add_option("--filename", filename, "File name prefix")->capture_default_str()->group("Output");

    SAMURAI_PARSE(argc, argv);

    if (!fs::exists(path))
    {
        fs::create_directory(path);
    }

    using MRMesh      = samurai::MRMesh<MRConfig>;
    using mrmesh_id_t = typename MRMesh::mesh_id_t;

    using UConfig = samurai::UniformConfig<dim>;
    using UMesh   = samurai::UniformMesh<UConfig>;

    const samurai::Box<double, dim> box({-1, -1, -1}, {1, 1, 1});
    MRMesh mrmesh{box, min_level, max_level, 0, 1};
    UMesh umesh{box, max_level, 0, 1};
    auto u       = init(mrmesh, test_case);
    auto u_exact = init(umesh, test_case);

    auto MRadaptation = samurai::make_MRAdapt(u);
    auto mra_config   = samurai::mra_config().regularity(2);
    MRadaptation(mra_config);

    auto level_ = samurai::make_scalar_field<std::size_t>("level", mrmesh);
    samurai::for_each_cell(mrmesh[mrmesh_id_t::cells],
                           [&](const auto& cell)
                           {
                               level_[cell] = cell.level;
                           });
    samurai::save(path, filename, mrmesh, u, level_);

    auto t1            = std::chrono::high_resolution_clock::now();
    auto u_reconstruct = reconstruction(u);
    auto t2            = std::chrono::high_resolution_clock::now();
    std::cout << "execution time " << std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count() << std::endl;

    auto error = samurai::make_scalar_field<double>("error", u_reconstruct.mesh());
    samurai::for_each_interval(u_reconstruct.mesh(),
                               [&](std::size_t level, const auto& i, const auto& index)
                               {
                                   auto j                = index[0];
                                   auto k                = index[1];
                                   error(level, i, j, k) = abs(u_reconstruct(level, i, j, k) - u_exact(level, i, j, k));
                               });
    samurai::save(path, fmt::format("uniform_{}", filename), u_reconstruct.mesh(), u_reconstruct, error);

    samurai::finalize();
    return 0;
}
