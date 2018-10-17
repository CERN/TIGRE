#!/bin/bash

rm -r build dist tigre.egg-info .idea 
rm _Ax.so _Atb.so _tvdenoising.so
rm tigre/Source/_Ax.cpp tigre/Source/_Atb.cpp
pip uninstall tigre