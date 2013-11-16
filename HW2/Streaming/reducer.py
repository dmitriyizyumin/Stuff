#!/usr/bin/env python

import sys

my_bin = None
my_count = 0

for line in sys.stdin:
    line = line.strip()
    this_bin,count=line.split("\t",1)
    count = int(count)

    if my_bin == this_bin:
        my_count += count
    else:
         if my_bin:
            (bin1,bin2,bin3,bin4) = my_bin.split(",",3)
            new_row = [bin1,bin2,bin3,bin4,my_count]
            print '%s,%s,%s,%s,%s' % (bin1,bin2,bin3,bin4,my_count)

         my_count = count
         my_bin = this_bin

if my_bin == this_bin:
    (bin1,bin2,bin3,bin4) = my_bin.split(",",3)
    new_row = [bin1,bin2,bin3,bin4,my_count]
    print '%s,%s,%s,%s,%s' % (bin1,bin2,bin3,bin4,my_count)
