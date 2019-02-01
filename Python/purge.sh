#!/bin/bash

rm -r build dist tigre.egg-info .idea 
rm _Ax.so _Atb.so _tvdenoising.so _minTV.so
rm /home/reuben/TIGRE/Python/tigre/Source/_AwminTV.cpp
rm tigre/Source/_Ax.cpp tigre/Source/_Atb.cpp tigre/Source/_tvdenoising.cpp tigre/Source/_minTV.cpp tigre/Source/_AwminTV.cpp
pip uninstall tigre
rm /home/reuben/TIGRE/Python/_AwminTV.so