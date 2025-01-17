# test-mdspan

## Different implementation

My first version is to index a one dimension array.

```cpp
class field3_1dp
{
protected:
    unsigned int Nx, Ny, Nz;

public:
    double* value = nullptr;
    field3_1dp(unsigned int _Nx, unsigned int _Ny, unsigned int _Nz)
        : Nx(_Nx)
        , Ny(_Ny)
        , Nz(_Nz)
    {
        value = new double[Nx * Ny * Nz];

        for (unsigned int i = 0; i < Nx * Ny * Nz; i++)
            value[i] = 0.;
    }

    ~field3_1dp() { delete[] value; }

    double& operator()(unsigned int i, unsigned int j, unsigned int k) { return value[i * Ny * Nz + j * Nz + k]; }

    int SizeX() { return Nx; }
    int SizeY() { return Ny; }
    int SizeZ() { return Nz; }
};
```

To accerlate the indexing to one dimension array, I use template to hard code the dimension numbers.

```cpp
template<unsigned int Nx, unsigned int Ny, unsigned int Nz>
class field3_1dp_t
{
public:
    double* value = nullptr;
    field3_1dp_t()
    {
        value = new double[Nx * Ny * Nz];

        for (unsigned int i = 0; i < Nx * Ny * Nz; i++)
            value[i] = 0.;
    }

    ~field3_1dp_t() { delete[] value; }

    double& operator()(unsigned int i, unsigned int j, unsigned int k) { return value[i * Ny * Nz + j * Nz + k]; }
};
```

But you can also store three dimension array in a class.

```cpp
class field3_3dp
{
protected:
    unsigned int Nx, Ny, Nz;

public:
    double*** value = nullptr;
    field3_3dp(unsigned int _Nx, unsigned int _Ny, unsigned int _Nz)
        : Nx(_Nx)
        , Ny(_Ny)
        , Nz(_Nz)
    {
        value = allocate_3d_array(Nx, Ny, Nz);
    }

    ~field3_3dp() { deallocate_3d_array(value, Nx, Ny); }

    double& operator()(unsigned int i, unsigned int j, unsigned int k) { return value[i][j][k]; }

    int SizeX() { return Nx; }
    int SizeY() { return Ny; }
    int SizeZ() { return Nz; }
};
```

However, three dimension array will lead to dicontigous memory.

So why not map a one dimension array to three dimension pointer?

```cpp
class field3_map
{
protected:
    unsigned int Nx, Ny, Nz;
    double*      value;
    double***    ptr3d;

public:
    field3_map(unsigned int _Nx, unsigned int _Ny, unsigned int _Nz)
        : Nx(_Nx)
        , Ny(_Ny)
        , Nz(_Nz)
        , value(new double[Nx * Ny * Nz])
    {
        ptr3d = new double**[Nx];
        for (unsigned int i = 0; i < Nx; ++i)
        {
            ptr3d[i] = new double*[Ny];
            for (unsigned int j = 0; j < Ny; ++j)
            {
                ptr3d[i][j] = value + (i * Ny + j) * Nz;
            }
        }
    }

    ~field3_map()
    {
        for (unsigned int i = 0; i < Nx; ++i)
        {
            delete[] ptr3d[i];
        }
        delete[] ptr3d;
        delete[] value;
    }

    double** operator[](unsigned int i) { return ptr3d[i]; }

    int SizeX() { return Nx; }
    int SizeY() { return Ny; }
    int SizeZ() { return Nz; }
};
```

Kokkos has `mdspan` to do this.

```cpp
class field3_mdspan
{
protected:
    unsigned int                                     Nx, Ny, Nz;
    double*                                          value;
    Kokkos::mdspan<double, Kokkos::dextents<int, 3>> m_mdspan;

public:
    field3_mdspan(unsigned int _Nx, unsigned int _Ny, unsigned int _Nz)
        : Nx(_Nx)
        , Ny(_Ny)
        , Nz(_Nz)
        , value(new double[Nx * Ny * Nz])
        , m_mdspan(value, Nx, Ny, Nz)
    {}

    ~field3_mdspan() { delete[] value; }

    double& operator()(int i, int j, int k) { return m_mdspan[i, j, k]; }

    int SizeX() { return Nx; }
    int SizeY() { return Ny; }
    int SizeZ() { return Nz; }
};
```

Comparsion is conducted on convection funciton, which is widely used in CFD.

Each class have its own funciton due to parameters' type are different. One example is:

```cpp
// Function to calculate convection term u * ∇u + v * ∇v + w * ∇w
void calculate_convection_3dp_raw(double***    u,
                                  double***    v,
                                  double***    w,
                                  double***    conv_u,
                                  double***    conv_v,
                                  double***    conv_w,
                                  unsigned int Nx,
                                  unsigned int Ny,
                                  unsigned int Nz)
{
    for (unsigned int i = 1; i < Nx - 1; ++i)
    {
        for (unsigned int j = 1; j < Ny - 1; ++j)
        {
            for (unsigned int k = 1; k < Nz - 1; ++k)
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
```

Measurment is:

```cpp
// Measure the time for 100 iterations of convection calculation
const int num_iterations = 20;
auto      start_time     = std::chrono::high_resolution_clock::now();

for (int iter = 0; iter < num_iterations; ++iter)
{
    calculate_convection_1dp(u, v, w, conv_u, conv_v, conv_w);
    std::cout << "iter = " << iter << std::endl;
}

auto                          end_time     = std::chrono::high_resolution_clock::now();
std::chrono::duration<double> elapsed_time = end_time - start_time;
```

Here is comparsion result.

## Data Scale

Compiler: MSVC 19.39.33520.0

<table>
    <tr>
        <td align="center">class</td>
        <td colspan="6" align="center">time(s)</td>
    </tr>
    <tr>
        <td align="center">scale</td>
        <td align="center">field3_1dp</td>
        <td align="center">field3_1dp_t</td>
        <td align="center">field3_3dp</td>
        <td align="center">field3_r</td>
        <td align="center">field3_map</td>
        <td align="center">field3_mdspan</td>
    </tr>
    <tr>
        <td align="center">134217728</td>
        <td align="center">5.57658</td>
        <td align="center">3.96567</td>
        <td align="center">1.19348</td>
        <td align="center">1.0934</td>
        <td align="center">4.90802</td>
        <td align="center">2.35261</td>
    </tr>
    <tr>
        <td align="center">89915392</td>
        <td align="center">1.68941</td>
        <td align="center">1.18209</td>
        <td align="center">0.609244</td>
        <td align="center">0.620306</td>
        <td align="center">1.41848</td>
        <td align="center">0.780205</td>
    </tr>
    <tr>
        <td align="center">56623104</td>
        <td align="center">1.05672</td>
        <td align="center">0.761333</td>
        <td align="center">0.401177</td>
        <td align="center">0.408533</td>
        <td align="center">0.906887</td>
        <td align="center">0.480399</td>
    </tr>
    <tr>
        <td align="center">32768000</td>
        <td align="center">0.609183</td>
        <td align="center">0.456157</td>
        <td align="center">0.228897</td>
        <td align="center">0.228498</td>
        <td align="center">0.544337</td>
        <td align="center">0.286498</td>
    </tr>
    <tr>
        <td align="center">16777216</td>
        <td align="center">0.310148</td>
        <td align="center">0.225518</td>
        <td align="center">0.120197</td>
        <td align="center">0.122801</td>
        <td align="center">0.266692</td>
        <td align="center">0.155674</td>
    </tr>
    <tr>
        <td align="center">2097152</td>
        <td align="center">0.0377825</td>
        <td align="center">0.0300964</td>
        <td align="center">0.0156596</td>
        <td align="center">0.0155079</td>
        <td align="center">0.0317008</td>
        <td align="center">0.0169362</td>
    </tr>
    <tr>
        <td align="center">262144</td>
        <td align="center">0.00506641</td>
        <td align="center">0.00398439</td>
        <td align="center">0.00241431</td>
        <td align="center">0.00260458</td>
        <td align="center">0.00445506</td>
        <td align="center">0.00317016</td>
    </tr>
</table>

![](./docs/result.svg)

## Platform and Compiler

As the answer saying:

[https://www.reddit.com/r/cpp_questions/comments/1e8i8m1/performance_test_for_multiple_dimension_array_and/](https://www.reddit.com/r/cpp_questions/comments/1e8i8m1/performance_test_for_multiple_dimension_array_and/)

I replace all indexing and accessing type to `std::size_t`, add `-march=native` option. Then repeat the test.

```shell
g++ .\main.cpp -O3 -I3rdparty/mdspan/include -march=native -o main.exe
```

```shell
clang++ .\main.cpp -O3 -I3rdparty/mdspan/include -march=native -o main.exe
```

I can get the same result that all implementation are similar when using clang in Linux.

<table>
    <tr>
        <td rowspan="2" align="center">platform</td>
        <td rowspan="2" align="center">complier</td>
        <td colspan="6" align="center">time(s)</td>
    </tr>
    <tr>
        <td align="center">field3_1dp</td>
        <td align="center">field3_1dp_t</td>
        <td align="center">field3_3dp</td>
        <td align="center">field3_r</td>
        <td align="center">field3_map</td>
        <td align="center">field3_mdspan</td>
    </tr>
    <tr>
        <td rowspan="3" align="center">windows</td>
        <td align="center">msvc</td>
        <td align="center">6.1776</td>
        <td align="center">3.916</td>
        <td align="center">1.37877</td>
        <td align="center">1.03</td>
        <td align="center">6.43983</td>
        <td align="center">6.26672</td>
    </tr>
    <tr>
        <td align="center">g++</td>
        <td align="center">5.63062</td>
        <td align="center">4.18603</td>
        <td align="center">1.12186</td>
        <td align="center">1.06282</td>
        <td align="center">5.05098</td>
        <td align="center">5.11213</td>
    </tr>
    <tr>
        <td align="center">clang++</td>
        <td align="center">3.14699</td>
        <td align="center">1.68778</td>
        <td align="center">1.19851</td>
        <td align="center">1.13926</td>
        <td align="center">1.99766</td>
        <td align="center">1.93895</td>
    </tr>
    <tr>
        <td rowspan="2" align="center">linux</td>
        <td align="center">g++</td>
        <td align="center">1.90795</td>
        <td align="center">1.86067</td>
        <td align="center">1.42331</td>
        <td align="center">1.42438</td>
        <td align="center"> 1.87503</td>
        <td align="center">1.92114</td>
    </tr>
    <tr>
        <td align="center">clang++</td>
        <td align="center">1.03195</td>
        <td align="center">0.96762</td>
        <td align="center">1.10214</td>
        <td align="center">1.0981</td>
        <td align="center">1.04556</td>
        <td align="center">1.0522</td>
    </tr>
</table>