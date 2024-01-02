#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright © 2024 Ye Chang yech1990@gmail.com
# Distributed under terms of the GNU license.
#
# Created: 2024-01-01 17:07

import gzip
import rich_click as click


def open_file(filename, mode="r"):
    if filename.endswith(".gz"):
        if mode == "w":
            return gzip.open(filename, "wt")
        return gzip.open(filename, "rt")
    else:
        if mode == "w":
            return click.open_file(filename, "w")
        return click.open_file(filename, "r")
