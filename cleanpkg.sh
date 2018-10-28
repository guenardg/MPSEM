#!/bin/bash
rm -f *~
rm -f MPSEM/*~
rm -f MPSEM/man/*~
rm -f MPSEM/R/*~
rm -f MPSEM/src/*~
rm -f MPSEM/src/*.so
rm -f MPSEM/src/*.o
rm -f MPSEM/src/*.rds
cd MPSEM && find -type f \( -not -name "MD5" \) -exec md5sum '{}' \; > MD5
