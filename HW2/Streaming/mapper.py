#!/usr/bin/env python

import sys
import math

for line in sys.stdin:
    nums = line.split()
    x_lo = math.floor(float(nums[0])*10)/10
    x_hi = math.ceil(float(nums[0])*10)/10
    y_lo = math.floor(float(nums[1])*10)/10
    y_hi = math.ceil(float(nums[1])*10)/10
    print '%s,%s,%s,%s\t%s' % (x_lo,x_hi,y_lo,y_hi,1)

