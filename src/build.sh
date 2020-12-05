#!/bin/sh

swig -c++ -python -I/public/home/yinghuang/software/anaconda3/include/python3.7m -I/public/home/yinghuang/software/anaconda3/lib/python3.7/site-packages/mpi4py/include  mpiplus.i

mpic++ -I/public/home/yinghuang/software/anaconda3/include/python3.7m -I/public/home/yinghuang/software/anaconda3/lib/python3.7/site-packages/mpi4py/include  -o _mpiplus.so mpiplus.cpp mpiplus_wrap.cxx -fPIC -shared -lpthread -ldl -lutil
