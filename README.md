THE CAMERA TEMPORAL ANALYSIS DATA PROJECT

The source file in which you can modify/optimize the analysis is main.cpp

The file functions.h has to contain only functions used for the analysis program (called by main.cpp)
If the file functions.h exists, it should be added in the CMakeList.txt by the command

COMPILATION WITH CMAKE

$ cmake -H. -Bbuild
$ cmake --build build -- -j8

The project will be built in a folder named 'build' and you can execute the project with the command line:

./<PROJECT_NAME>

VECTORIZATION OF LQUA FILES FOR COMPUTING THE X Y Z COORDINATES

1- Install Docker Desktop for Mac

2- Loading the image
docker load -i vectorize-lqua.tar.gz (or vectorize-lqua.tar)

3- Run the script vectorize.lqua:

docker run --rm -ti -v pwd:/media vectorize-lqua vectorize-lqua Path_to_kqua Path_to_lqua_1 Path_to_lqua_2 ...

pwd correspond to the path of the mounting point of Docker. It will be the root of Docker ! The path to the lqua and kqua files have to be relative to the directory mounting point.

Advice : The mounting point of Docker should be the directory where .kqua and .lqua are stored in order to facilitate the autocompletion.
