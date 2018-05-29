#!/usr/bin/env python
# encoding: utf-8

import sys
sys.path.append('build/lib.macosx-10.7-x86_64-3.6/')
import linalg

v = linalg.Vector([1,2,3,4])
for i in range(100):
    v += v # print(v)
    m = linalg.Matrix([[1],[2],[3]])

print(m.getNRowVec(0))
