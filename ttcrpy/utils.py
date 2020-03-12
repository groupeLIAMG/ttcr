# -*- coding: utf-8 -*-

"""
    Created on Tue Jun 21 20:55:29 2016

    @author: giroux

    Copyright 2016 Bernard Giroux
    email: bernard.giroux@ete.inrs.ca

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it /will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""

import dis
import inspect
import sys


def nargout():
    """
    Return how many values the caller is expecting

    taken from
    http://stackoverflow.com/questions/16488872/python-check-for-the-number-of-output-arguments-a-function-is-called-with
    """
    if sys.version_info[1] > 5:
        off = 1
    else:
        off = 0
    f = inspect.currentframe()
    f = f.f_back.f_back
    c = f.f_code
    i = f.f_lasti
    bytecode = c.co_code
    instruction = bytecode[i+3-off]
    if instruction == dis.opmap['UNPACK_SEQUENCE']:
        howmany = bytecode[i+4-off]
        return howmany
    elif instruction == dis.opmap['POP_TOP']:
        return 0
    return 1
