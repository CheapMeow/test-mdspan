#include <chrono>
#include <cstdlib>
#include <iostream>
#include <mdspan/mdspan.hpp>

// Function to allocate a 3D array
double*** allocate_3d_array(std::size_t Nx, std::size_t Ny, std::size_t Nz)
{
    double*** array = new double**[Nx];
    for (std::size_t i = 0; i < Nx; ++i)
    {
        array[i] = new double*[Ny];
        for (std::size_t j = 0; j < Ny; ++j)
        {
            array[i][j] = new double[Nz];
        }
    }
    return array;
}

// Function to deallocate a 3D array
void deallocate_3d_array(double*** array, std::size_t Nx, std::size_t Ny)
{
    for (std::size_t i = 0; i < Nx; ++i)
    {
        for (std::size_t j = 0; j < Ny; ++j)
        {
            delete[] array[i][j];
        }
        delete[] array[i];
    }
    delete[] array;
}

class field3_1dp
{
protected:
    std::size_t Nx, Ny, Nz;

public:
    double* value = nullptr;
    field3_1dp(std::size_t _Nx, std::size_t _Ny, std::size_t _Nz)
        : Nx(_Nx)
        , Ny(_Ny)
        , Nz(_Nz)
    {
        value = new double[Nx * Ny * Nz];

        for (std::size_t i = 0; i < Nx * Ny * Nz; i++)
            value[i] = 0.;
    }

    ~field3_1dp() { delete[] value; }

    double& operator()(std::size_t i, std::size_t j, std::size_t k) { return value[i * Ny * Nz + j * Nz + k]; }

    std::size_t SizeX() { return Nx; }
    std::size_t SizeY() { return Ny; }
    std::size_t SizeZ() { return Nz; }
};

// Function to calculate convection term u * ∇u + v * ∇v + w * ∇w
void calculate_convection_1dp(field3_1dp& u,
                              field3_1dp& v,
                              field3_1dp& w,
                              field3_1dp& conv_u,
                              field3_1dp& conv_v,
                              field3_1dp& conv_w)
{
    for (std::size_t i = 1; i < u.SizeX() - 1; ++i)
    {
        for (std::size_t j = 1; j < u.SizeY() - 1; ++j)
        {
            for (std::size_t k = 1; k < u.SizeZ() - 1; ++k)
            {
                double dudx = (u(i + 1, j, k) - u(i - 1, j, k)) / 2.0;
                double dudy = (u(i, j + 1, k) - u(i, j - 1, k)) / 2.0;
                double dudz = (u(i, j, k + 1) - u(i, j, k - 1)) / 2.0;

                double dvdx = (v(i + 1, j, k) - v(i - 1, j, k)) / 2.0;
                double dvdy = (v(i, j + 1, k) - v(i, j - 1, k)) / 2.0;
                double dvdz = (v(i, j, k + 1) - v(i, j, k - 1)) / 2.0;

                double dwdx = (w(i + 1, j, k) - w(i - 1, j, k)) / 2.0;
                double dwdy = (w(i, j + 1, k) - w(i, j - 1, k)) / 2.0;
                double dwdz = (w(i, j, k + 1) - w(i, j, k - 1)) / 2.0;

                conv_u(i, j, k) = u(i, j, k) * dudx + v(i, j, k) * dudy + w(i, j, k) * dudz;
                conv_v(i, j, k) = u(i, j, k) * dvdx + v(i, j, k) * dvdy + w(i, j, k) * dvdz;
                conv_w(i, j, k) = u(i, j, k) * dwdx + v(i, j, k) * dwdy + w(i, j, k) * dwdz;
            }
        }
    }
}

template<std::size_t Nx, std::size_t Ny, std::size_t Nz>
class field3_1dp_t
{
public:
    double* value = nullptr;
    field3_1dp_t()
    {
        value = new double[Nx * Ny * Nz];

        for (std::size_t i = 0; i < Nx * Ny * Nz; i++)
            value[i] = 0.;
    }

    ~field3_1dp_t() { delete[] value; }

    double& operator()(std::size_t i, std::size_t j, std::size_t k) { return value[i * Ny * Nz + j * Nz + k]; }
};

template<std::size_t Nx, std::size_t Ny, std::size_t Nz>
void calculate_convection_1dp_t(field3_1dp_t<Nx, Ny, Nz>& u,
                                field3_1dp_t<Nx, Ny, Nz>& v,
                                field3_1dp_t<Nx, Ny, Nz>& w,
                                field3_1dp_t<Nx, Ny, Nz>& conv_u,
                                field3_1dp_t<Nx, Ny, Nz>& conv_v,
                                field3_1dp_t<Nx, Ny, Nz>& conv_w)
{
    for (std::size_t i = 1; i < Nx - 1; ++i)
    {
        for (std::size_t j = 1; j < Ny - 1; ++j)
        {
            for (std::size_t k = 1; k < Nz - 1; ++k)
            {
                double dudx = (u(i + 1, j, k) - u(i - 1, j, k)) / 2.0;
                double dudy = (u(i, j + 1, k) - u(i, j - 1, k)) / 2.0;
                double dudz = (u(i, j, k + 1) - u(i, j, k - 1)) / 2.0;

                double dvdx = (v(i + 1, j, k) - v(i - 1, j, k)) / 2.0;
                double dvdy = (v(i, j + 1, k) - v(i, j - 1, k)) / 2.0;
                double dvdz = (v(i, j, k + 1) - v(i, j, k - 1)) / 2.0;

                double dwdx = (w(i + 1, j, k) - w(i - 1, j, k)) / 2.0;
                double dwdy = (w(i, j + 1, k) - w(i, j - 1, k)) / 2.0;
                double dwdz = (w(i, j, k + 1) - w(i, j, k - 1)) / 2.0;

                conv_u(i, j, k) = u(i, j, k) * dudx + v(i, j, k) * dudy + w(i, j, k) * dudz;
                conv_v(i, j, k) = u(i, j, k) * dvdx + v(i, j, k) * dvdy + w(i, j, k) * dvdz;
                conv_w(i, j, k) = u(i, j, k) * dwdx + v(i, j, k) * dwdy + w(i, j, k) * dwdz;
            }
        }
    }
}

class field3_3dp
{
protected:
    std::size_t Nx, Ny, Nz;

public:
    double*** value = nullptr;
    field3_3dp(std::size_t _Nx, std::size_t _Ny, std::size_t _Nz)
        : Nx(_Nx)
        , Ny(_Ny)
        , Nz(_Nz)
    {
        value = allocate_3d_array(Nx, Ny, Nz);
    }

    ~field3_3dp() { deallocate_3d_array(value, Nx, Ny); }

    double& operator()(std::size_t i, std::size_t j, std::size_t k) { return value[i][j][k]; }

    std::size_t SizeX() { return Nx; }
    std::size_t SizeY() { return Ny; }
    std::size_t SizeZ() { return Nz; }
};

void calculate_convection_2(field3_3dp& u,
                            field3_3dp& v,
                            field3_3dp& w,
                            field3_3dp& conv_u,
                            field3_3dp& conv_v,
                            field3_3dp& conv_w)
{
    for (std::size_t i = 1; i < u.SizeX() - 1; ++i)
    {
        for (std::size_t j = 1; j < u.SizeY() - 1; ++j)
        {
            for (std::size_t k = 1; k < u.SizeZ() - 1; ++k)
            {
                double dudx = (u(i + 1, j, k) - u(i - 1, j, k)) / 2.0;
                double dudy = (u(i, j + 1, k) - u(i, j - 1, k)) / 2.0;
                double dudz = (u(i, j, k + 1) - u(i, j, k - 1)) / 2.0;

                double dvdx = (v(i + 1, j, k) - v(i - 1, j, k)) / 2.0;
                double dvdy = (v(i, j + 1, k) - v(i, j - 1, k)) / 2.0;
                double dvdz = (v(i, j, k + 1) - v(i, j, k - 1)) / 2.0;

                double dwdx = (w(i + 1, j, k) - w(i - 1, j, k)) / 2.0;
                double dwdy = (w(i, j + 1, k) - w(i, j - 1, k)) / 2.0;
                double dwdz = (w(i, j, k + 1) - w(i, j, k - 1)) / 2.0;

                conv_u(i, j, k) = u(i, j, k) * dudx + v(i, j, k) * dudy + w(i, j, k) * dudz;
                conv_v(i, j, k) = u(i, j, k) * dvdx + v(i, j, k) * dvdy + w(i, j, k) * dvdz;
                conv_w(i, j, k) = u(i, j, k) * dwdx + v(i, j, k) * dwdy + w(i, j, k) * dwdz;
            }
        }
    }
}

class field3_map
{
protected:
    std::size_t Nx, Ny, Nz;
    double*     value;
    double***   ptr3d;

public:
    field3_map(std::size_t _Nx, std::size_t _Ny, std::size_t _Nz)
        : Nx(_Nx)
        , Ny(_Ny)
        , Nz(_Nz)
        , value(new double[Nx * Ny * Nz])
    {
        ptr3d = new double**[Nx];
        for (std::size_t i = 0; i < Nx; ++i)
        {
            ptr3d[i] = new double*[Ny];
            for (std::size_t j = 0; j < Ny; ++j)
            {
                ptr3d[i][j] = value + (i * Ny + j) * Nz;
            }
        }
    }

    ~field3_map()
    {
        for (std::size_t i = 0; i < Nx; ++i)
        {
            delete[] ptr3d[i];
        }
        delete[] ptr3d;
        delete[] value;
    }

    double** operator[](std::size_t i) { return ptr3d[i]; }

    std::size_t SizeX() { return Nx; }
    std::size_t SizeY() { return Ny; }
    std::size_t SizeZ() { return Nz; }
};

void calculate_convection_map(field3_map& u,
                              field3_map& v,
                              field3_map& w,
                              field3_map& conv_u,
                              field3_map& conv_v,
                              field3_map& conv_w)
{
    for (std::size_t i = 1; i < u.SizeX() - 1; ++i)
    {
        for (std::size_t j = 1; j < u.SizeY() - 1; ++j)
        {
            for (std::size_t k = 1; k < u.SizeZ() - 1; ++k)
            {
                double dudx = (u[i + 1][j][k] - u[i - 1][j][k]) / 2.0;
                double dudy = (u[i][j + 1][k] - u[i][j - 1][k]) / 2.0;
                double dudz = (u[i][j][k + 1] - u[i][j][k - 1]) / 2.0;

                double dvdx = (v[i + 1][j][k] - v[i - 1][j][k]) / 2.0;
                double dvdy = (v[i][j + 1][k] - v[i][j - 1][k]) / 2.0;
                double dvdz = (v[i][j][k + 1] - v[i][j][k - 1]) / 2.0;

                double dwdx = (w[i + 1][j][k] - w[i - 1][j][k]) / 2.0;
                double dwdy = (w[i][j + 1][k] - w[i][j - 1][k]) / 2.0;
                double dwdz = (w[i][j][k + 1] - w[i][j][k - 1]) / 2.0;

                conv_u[i][j][k] = u[i][j][k] * dudx + v[i][j][k] * dudy + w[i][j][k] * dudz;
                conv_v[i][j][k] = u[i][j][k] * dvdx + v[i][j][k] * dvdy + w[i][j][k] * dvdz;
                conv_w[i][j][k] = u[i][j][k] * dwdx + v[i][j][k] * dwdy + w[i][j][k] * dwdz;
            }
        }
    }
}

class field3_mdspan
{
protected:
    std::size_t                                              Nx, Ny, Nz;
    double*                                                  value;
    Kokkos::mdspan<double, Kokkos::dextents<std::size_t, 3>> m_mdspan;

public:
    field3_mdspan(std::size_t _Nx, std::size_t _Ny, std::size_t _Nz)
        : Nx(_Nx)
        , Ny(_Ny)
        , Nz(_Nz)
        , value(new double[Nx * Ny * Nz])
        , m_mdspan(value, Nx, Ny, Nz)
    {}

    ~field3_mdspan() { delete[] value; }

    double& operator()(std::size_t i, std::size_t j, std::size_t k) { return m_mdspan(i, j, k); }

    std::size_t SizeX() { return Nx; }
    std::size_t SizeY() { return Ny; }
    std::size_t SizeZ() { return Nz; }
};

void calculate_convection_mdspan(field3_mdspan& u,
                                 field3_mdspan& v,
                                 field3_mdspan& w,
                                 field3_mdspan& conv_u,
                                 field3_mdspan& conv_v,
                                 field3_mdspan& conv_w)
{
    for (std::size_t i = 1; i < u.SizeX() - 1; ++i)
    {
        for (std::size_t j = 1; j < u.SizeY() - 1; ++j)
        {
            for (std::size_t k = 1; k < u.SizeZ() - 1; ++k)
            {
                double dudx = (u(i + 1, j, k) - u(i - 1, j, k)) / 2.0;
                double dudy = (u(i, j + 1, k) - u(i, j - 1, k)) / 2.0;
                double dudz = (u(i, j, k + 1) - u(i, j, k - 1)) / 2.0;

                double dvdx = (v(i + 1, j, k) - v(i - 1, j, k)) / 2.0;
                double dvdy = (v(i, j + 1, k) - v(i, j - 1, k)) / 2.0;
                double dvdz = (v(i, j, k + 1) - v(i, j, k - 1)) / 2.0;

                double dwdx = (w(i + 1, j, k) - w(i - 1, j, k)) / 2.0;
                double dwdy = (w(i, j + 1, k) - w(i, j - 1, k)) / 2.0;
                double dwdz = (w(i, j, k + 1) - w(i, j, k - 1)) / 2.0;

                conv_u(i, j, k) = u(i, j, k) * dudx + v(i, j, k) * dudy + w(i, j, k) * dudz;
                conv_v(i, j, k) = u(i, j, k) * dvdx + v(i, j, k) * dvdy + w(i, j, k) * dvdz;
                conv_w(i, j, k) = u(i, j, k) * dwdx + v(i, j, k) * dwdy + w(i, j, k) * dwdz;
            }
        }
    }
}

// Function to calculate convection term u * ∇u + v * ∇v + w * ∇w
void calculate_convection_3dp_raw(double***   u,
                                  double***   v,
                                  double***   w,
                                  double***   conv_u,
                                  double***   conv_v,
                                  double***   conv_w,
                                  std::size_t Nx,
                                  std::size_t Ny,
                                  std::size_t Nz)
{
    for (std::size_t i = 1; i < Nx - 1; ++i)
    {
        for (std::size_t j = 1; j < Ny - 1; ++j)
        {
            for (std::size_t k = 1; k < Nz - 1; ++k)
            {
                double dudx = (u[i + 1][j][k] - u[i - 1][j][k]) / 2.0;
                double dudy = (u[i][j + 1][k] - u[i][j - 1][k]) / 2.0;
                double dudz = (u[i][j][k + 1] - u[i][j][k - 1]) / 2.0;

                double dvdx = (v[i + 1][j][k] - v[i - 1][j][k]) / 2.0;
                double dvdy = (v[i][j + 1][k] - v[i][j - 1][k]) / 2.0;
                double dvdz = (v[i][j][k + 1] - v[i][j][k - 1]) / 2.0;

                double dwdx = (w[i + 1][j][k] - w[i - 1][j][k]) / 2.0;
                double dwdy = (w[i][j + 1][k] - w[i][j - 1][k]) / 2.0;
                double dwdz = (w[i][j][k + 1] - w[i][j][k - 1]) / 2.0;

                conv_u[i][j][k] = u[i][j][k] * dudx + v[i][j][k] * dudy + w[i][j][k] * dudz;
                conv_v[i][j][k] = u[i][j][k] * dvdx + v[i][j][k] * dvdy + w[i][j][k] * dvdz;
                conv_w[i][j][k] = u[i][j][k] * dwdx + v[i][j][k] * dwdy + w[i][j][k] * dwdz;
            }
        }
    }
}

int main()
{
    const std::size_t Nx = 512, Ny = 512, Nz = 512;

    {
        field3_1dp u(Nx, Ny, Nz);
        field3_1dp v(Nx, Ny, Nz);
        field3_1dp w(Nx, Ny, Nz);
        field3_1dp conv_u(Nx, Ny, Nz);
        field3_1dp conv_v(Nx, Ny, Nz);
        field3_1dp conv_w(Nx, Ny, Nz);

        // Assign some values to u, v, and w for testing
        for (std::size_t i = 0; i < Nx; ++i)
        {
            for (std::size_t j = 0; j < Ny; ++j)
            {
                for (std::size_t k = 0; k < Nz; ++k)
                {
                    u(i, j, k) = static_cast<double>(i + j + k);
                    v(i, j, k) = static_cast<double>(i + j + k);
                    w(i, j, k) = static_cast<double>(i + j + k);
                }
            }
        }

        // Measure the time for 100 iterations of convection calculation
        const std::size_t num_iterations = 20;
        auto              start_time     = std::chrono::high_resolution_clock::now();

        for (std::size_t iter = 0; iter < num_iterations; ++iter)
        {
            calculate_convection_1dp(u, v, w, conv_u, conv_v, conv_w);
            std::cout << "iter = " << iter << std::endl;
        }

        auto                          end_time     = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed_time = end_time - start_time;

        std::cout << "Average time per iteration for field3_1dp: " << (elapsed_time.count() / num_iterations)
                  << " seconds" << std::endl;
    }

    {
        field3_1dp_t<Nx, Ny, Nz> u;
        field3_1dp_t<Nx, Ny, Nz> v;
        field3_1dp_t<Nx, Ny, Nz> w;
        field3_1dp_t<Nx, Ny, Nz> conv_u;
        field3_1dp_t<Nx, Ny, Nz> conv_v;
        field3_1dp_t<Nx, Ny, Nz> conv_w;

        // Assign some values to u, v, and w for testing
        for (std::size_t i = 0; i < Nx; ++i)
        {
            for (std::size_t j = 0; j < Ny; ++j)
            {
                for (std::size_t k = 0; k < Nz; ++k)
                {
                    u(i, j, k) = static_cast<double>(i + j + k);
                    v(i, j, k) = static_cast<double>(i + j + k);
                    w(i, j, k) = static_cast<double>(i + j + k);
                }
            }
        }

        // Measure the time for 100 iterations of convection calculation
        const std::size_t num_iterations = 20;
        auto              start_time     = std::chrono::high_resolution_clock::now();

        for (std::size_t iter = 0; iter < num_iterations; ++iter)
        {
            calculate_convection_1dp_t(u, v, w, conv_u, conv_v, conv_w);
            std::cout << "iter = " << iter << std::endl;
        }

        auto                          end_time     = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed_time = end_time - start_time;

        std::cout << "Average time per iteration for field3_1dp_t: " << (elapsed_time.count() / num_iterations)
                  << " seconds" << std::endl;
    }

    {
        field3_3dp u(Nx, Ny, Nz);
        field3_3dp v(Nx, Ny, Nz);
        field3_3dp w(Nx, Ny, Nz);
        field3_3dp conv_u(Nx, Ny, Nz);
        field3_3dp conv_v(Nx, Ny, Nz);
        field3_3dp conv_w(Nx, Ny, Nz);

        // Assign some values to u, v, and w for testing
        for (std::size_t i = 0; i < Nx; ++i)
        {
            for (std::size_t j = 0; j < Ny; ++j)
            {
                for (std::size_t k = 0; k < Nz; ++k)
                {
                    u(i, j, k) = static_cast<double>(i + j + k);
                    v(i, j, k) = static_cast<double>(i + j + k);
                    w(i, j, k) = static_cast<double>(i + j + k);
                }
            }
        }

        // Measure the time for 100 iterations of convection calculation
        const std::size_t num_iterations = 20;
        auto              start_time     = std::chrono::high_resolution_clock::now();

        for (std::size_t iter = 0; iter < num_iterations; ++iter)
        {
            calculate_convection_2(u, v, w, conv_u, conv_v, conv_w);
            std::cout << "iter = " << iter << std::endl;
        }

        auto                          end_time     = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed_time = end_time - start_time;

        std::cout << "Average time per iteration for field3_3dp: " << (elapsed_time.count() / num_iterations)
                  << " seconds" << std::endl;
    }

    {
        double*** u      = allocate_3d_array(Nx, Ny, Nz);
        double*** v      = allocate_3d_array(Nx, Ny, Nz);
        double*** w      = allocate_3d_array(Nx, Ny, Nz);
        double*** conv_u = allocate_3d_array(Nx, Ny, Nz);
        double*** conv_v = allocate_3d_array(Nx, Ny, Nz);
        double*** conv_w = allocate_3d_array(Nx, Ny, Nz);

        // Assign some values to u, v, and w for testing
        for (std::size_t i = 0; i < Nx; ++i)
        {
            for (std::size_t j = 0; j < Ny; ++j)
            {
                for (std::size_t k = 0; k < Nz; ++k)
                {
                    u[i][j][k] = static_cast<double>(i + j + k);
                    v[i][j][k] = static_cast<double>(i + j + k);
                    w[i][j][k] = static_cast<double>(i + j + k);
                }
            }
        }

        // Measure the time for 100 iterations of convection calculation
        const std::size_t num_iterations = 20;
        auto              start_time     = std::chrono::high_resolution_clock::now();

        for (std::size_t iter = 0; iter < num_iterations; ++iter)
        {
            calculate_convection_3dp_raw(u, v, w, conv_u, conv_v, conv_w, Nx, Ny, Nz);
            std::cout << "iter = " << iter << std::endl;
        }

        auto                          end_time     = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed_time = end_time - start_time;

        std::cout << "Average time per iteration for field3_r: " << (elapsed_time.count() / num_iterations)
                  << " seconds" << std::endl;

        // Deallocate the 3D arrays
        deallocate_3d_array(u, Nx, Ny);
        deallocate_3d_array(v, Nx, Ny);
        deallocate_3d_array(w, Nx, Ny);
        deallocate_3d_array(conv_u, Nx, Ny);
        deallocate_3d_array(conv_v, Nx, Ny);
        deallocate_3d_array(conv_w, Nx, Ny);
    }

    {
        field3_map u(Nx, Ny, Nz);
        field3_map v(Nx, Ny, Nz);
        field3_map w(Nx, Ny, Nz);
        field3_map conv_u(Nx, Ny, Nz);
        field3_map conv_v(Nx, Ny, Nz);
        field3_map conv_w(Nx, Ny, Nz);

        // Assign some values to u, v, and w for testing
        for (std::size_t i = 0; i < Nx; ++i)
        {
            for (std::size_t j = 0; j < Ny; ++j)
            {
                for (std::size_t k = 0; k < Nz; ++k)
                {
                    u[i][j][k] = static_cast<double>(i + j + k);
                    v[i][j][k] = static_cast<double>(i + j + k);
                    w[i][j][k] = static_cast<double>(i + j + k);
                }
            }
        }

        // Measure the time for 100 iterations of convection calculation
        const std::size_t num_iterations = 20;
        auto              start_time     = std::chrono::high_resolution_clock::now();

        for (std::size_t iter = 0; iter < num_iterations; ++iter)
        {
            calculate_convection_map(u, v, w, conv_u, conv_v, conv_w);
            std::cout << "iter = " << iter << std::endl;
        }

        auto                          end_time     = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed_time = end_time - start_time;

        std::cout << "Average time per iteration for field3_map: " << (elapsed_time.count() / num_iterations)
                  << " seconds" << std::endl;
    }

    {
        field3_mdspan u(Nx, Ny, Nz);
        field3_mdspan v(Nx, Ny, Nz);
        field3_mdspan w(Nx, Ny, Nz);
        field3_mdspan conv_u(Nx, Ny, Nz);
        field3_mdspan conv_v(Nx, Ny, Nz);
        field3_mdspan conv_w(Nx, Ny, Nz);

        // Assign some values to u, v, and w for testing
        for (std::size_t i = 0; i < Nx; ++i)
        {
            for (std::size_t j = 0; j < Ny; ++j)
            {
                for (std::size_t k = 0; k < Nz; ++k)
                {
                    u(i, j, k) = static_cast<double>(i + j + k);
                    v(i, j, k) = static_cast<double>(i + j + k);
                    w(i, j, k) = static_cast<double>(i + j + k);
                }
            }
        }

        // Measure the time for 100 iterations of convection calculation
        const std::size_t num_iterations = 20;
        auto              start_time     = std::chrono::high_resolution_clock::now();

        for (std::size_t iter = 0; iter < num_iterations; ++iter)
        {
            calculate_convection_mdspan(u, v, w, conv_u, conv_v, conv_w);
            std::cout << "iter = " << iter << std::endl;
        }

        auto                          end_time     = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed_time = end_time - start_time;

        std::cout << "Average time per iteration for field3_mdspan: " << (elapsed_time.count() / num_iterations)
                  << " seconds" << std::endl;
    }

    return 0;
}
