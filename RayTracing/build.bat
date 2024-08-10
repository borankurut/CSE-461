mkdir build
cd build
cmake .. 
cmake --build . --config Release
cd ..
build\Release\main.exe > image.ppm
.\image.ppm

