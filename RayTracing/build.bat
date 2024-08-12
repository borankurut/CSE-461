mkdir build
cd build
cmake -G "Ninja" -DCMAKE_BUILD_TYPE=Release ..
cmake --build .
cd ..
build\main.exe "scene_cat.xml" > scene_cat.ppm
build\main.exe "scene_chest.xml" > scene_chest.ppm
build\main.exe "scene_spheres.xml" > scene_spheres.ppm

