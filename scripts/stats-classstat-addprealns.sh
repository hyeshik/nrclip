#!/bin/sh
#
# Copyright (c) 2011 Seoul National University RNA Biology Laboratory.
# 
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
# 
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.
#
# - Hyeshik Chang <hyeshik@snu.ac.kr>
#

(cat $1;
 zcat $2 | \
  gawk -F'	' '
  /^[0-9]/ {
    if ($3 != "*" && previd != $1) {
      split($1, rn, "-");
      cls = $3;
      totalcnt[cls] += rn[2];
      previd = $1;
    }
  }
  END {
    for (cls in totalcnt) {
      printf "%s,%d\r\n", cls, totalcnt[cls];
    }
  }') | sort -t, -k2,2nr | \
gawk -F, '
BEGIN {
  print "class,reads,proportion\r";
}
/[0-9]/ {
  i += 1;
  label[i] = $1;
  readcnt[i] = $2;
  totalcnt += $2;
}
END {
  for (j = 1; j <= i; j++) {
    printf "%s,%d,%.10f\r\n", label[j], readcnt[j], readcnt[j] / totalcnt;
  }
}'
