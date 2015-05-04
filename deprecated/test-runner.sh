#!/bin/bash

set -e
./$1 > $1.out 2>&1
diff -q $1.out $1.cmp
