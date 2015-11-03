#!/bin/sh

progName="treeMatcherMacros";

g++ $progName.C $(root-config --cflags --libs) -Wall -O2 -o $progName.exe || exit 1

./$progName.exe || exit 1
