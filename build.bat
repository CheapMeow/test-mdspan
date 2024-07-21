REM Find MinGW-w64 gcc and g++
set GCC_PATH=E:\software\mingw64\bin\gcc.exe
set GXX_PATH=E:\software\mingw64\bin\g++.exe

REM Check if gcc exists
if not exist "%GCC_PATH%" (
    echo gcc not found at: %GCC_PATH%
    exit /b 1
) else (
    echo gcc found at: %GCC_PATH%
)

REM Check if g++ exists
if not exist "%GXX_PATH%" (
    echo g++ not found at: %GXX_PATH%
    exit /b 1
) else (
    echo g++ found at: %GXX_PATH%
)

REM Set CC and CXX environment variables
set CC=%GCC_PATH%
set CXX=%GXX_PATH%

cd build
cmake -G "MinGW Makefiles" -DCMAKE_BUILD_TYPE=Release -DCMAKE_C_COMPILER=%CC% -DCMAKE_CXX_COMPILER=%CXX% ..
mingw32-make -j8