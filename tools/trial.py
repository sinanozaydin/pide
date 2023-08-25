#!/usr/bin/env python3

import os,sys

core_path_ext = os.path.join(os.path.dirname(os.path.abspath(__file__)) , '../SEEL')
print(core_path_ext)

sys.path.append(core_path_ext)

import SEEL

