Tested on WSL Ubuntu 24.04.

- install libs required for geant

    > sudo apt-get install dpkg-dev cmake g++ gcc binutils libx11-dev libxpm-dev libxft-dev libxext-dev libxmu-dev python3 libxerces-c-dev qtbase5-dev

- clone https://github.com/aovsiannikova/327 

- build geant

    > wget https://geant4-data.web.cern.ch/releases/geant4.10.07.p02.tar.gz
    > tar zxf geant4.10.07.p02.tar.gz # will create geant4.10.07.p02 directory
    > mkdir geant-build
    > mkdir geant-install
    > cd geant-build
    > cmake -DCMAKE_INSTALL_PREFIX=../geant-install -DGEANT4_INSTALL_DATA=ON -DGEANT4_USE_OPENGL_X11=ON -DGEANT4_BUILD_MULTITHREADED=ON -DGEANT4_USE_QT=ON ../geant4.10.07.p02
    > patch -p1 < 327/geant_10.7.2.patch in the geant4.10.07.p02 directory.

    > make # or make -j8 to build faster in 8 threads
    > make install

- change cmake script `geant-install/lib/Geant4-10.7.2/Geant4PackageCache.cmake`

  Replace:

      geant4_set_and_check_package_variable(EXPAT_LIBRARY ""  "")

  with

      geant4_set_and_check_package_variable(EXPAT_LIBRARY "/usr/lib/x86_64-linux-gnu/libexpat.so.1" PATH "path to expat lib")

- build app:

    > cd 327
    > mkdir app-build
    > mkdir app-install
    > cd app-build
    > cmake -DCMAKE_INSTALL_PREFIX=../app-install -DGeant4_DIR=/path/to/geant-install/lib/Geant4-10.7.2 -DGEANT4_INSTALL_DATA=ON -DGEANT4_USE_OPENGL_X11=ON -DGEANT4_BUILD_MULTITHREADED=ON -DGEANT4_USE_QT=ON ..
    > make
    > make install


- run app:

    > source /path/to/geant-install/bin/geant4.sh
    > ../app-install/OpNovice2 <gamma.mac>



