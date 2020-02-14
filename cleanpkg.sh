#!/bin/bash
PKG='MPSEM'
rm -f *~
rm -f .*~
rm $PKG.Rcheck.tar.gz
tar cvzf $PKG.Rcheck.tar.gz $PKG.Rcheck
rm -rf $PKG.Rcheck
rm -f $PKG/*~
rm -f $PKG/man/*~
rm -f $PKG/R/*~
rm -f $PKG/src/*~
rm -f $PKG/src/*.so
rm -f $PKG/src/*.o
rm -f $PKG/src/*.rds
cd $PKG && find -type f \( -not -name "MD5" \) -exec md5sum '{}' \; > MD5
