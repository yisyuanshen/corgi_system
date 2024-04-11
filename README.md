# corgi_system

1. (Already done) Comment ```find_package(mip REQUIRED NO_MODULE)``` and ```find_package(mip REQUIRED NO_MODULE)``` in *CMakeLists.txt*

2. (Already done) Comment ```add_subdirectory(...)``` except for ```cpg, kinematic, fsm``` in *src/CMakeLists.txt*

3. Uncomment ```// import std.proto``` in all the proto files at *corgi_system/protos*

4. Build with the following command

>>```
>>cmake .. -DBUILD_WITH_ROS=OFF -DLOCAL_PACKAGE_PATH=$HOME/corgi_ws/build
>>```
