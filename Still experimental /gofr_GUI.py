#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 22 14:20:29 2023

@author: k
"""

import gofr_umd_or
import sys

def main(argv):
    print(argv)
    gofr_umd_or.main(argv)

if __name__=="__main__":
    main(sys.argv[1:])
