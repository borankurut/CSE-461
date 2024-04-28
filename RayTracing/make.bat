mkdir build
cd build
cmake .. 
cmake --build . --config Debug
cd ..
build\Debug\main.exe > image.ppm
.\image.ppm

