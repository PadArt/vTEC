#!/usr/bin/env python
# coding=utf8
#
# Copyright 2024 Artem Padokhin <padokhin@physics.msu.ru>
#
# This file is part of vTEC.
#
# vTEC is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# vTEC is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with vTEC.  If not, see <http://www.gnu.org/licenses/>.


import numpy as np
from datetime import datetime


def sec_of_day(time):
    day_start = time.replace(hour=0, minute=0, second=0, microsecond=0)
    return (time - day_start).total_seconds()
sec_of_day = np.vectorize(sec_of_day)


def sec_of_interval(time, time0):
    return (time - time0).total_seconds()
sec_of_interval = np.vectorize(sec_of_interval, excluded='time0')

