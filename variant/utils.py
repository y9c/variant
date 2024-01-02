#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright © 2024 Ye Chang yech1990@gmail.com
# Distributed under terms of the GNU license.
#
# Created: 2024-01-01 17:07

import gzip
import logging

import rich_click as click


def get_logger(name: str) -> logging.Logger:
    """global logging."""
    logger: logging.Logger = logging.getLogger(name)
    if not logger.handlers:
        handler: logging.StreamHandler = logging.StreamHandler()
        formatter: logging.Formatter = logging.Formatter(
            "%(asctime)s %(name)-12s %(levelname)-8s %(message)s"
        )
        handler.setFormatter(formatter)
        logger.addHandler(handler)
        #  logger.setLevel(logging.DEBUG)
        logger.setLevel(logging.INFO)
    return logger


LOGGER: logging.Logger = get_logger(__name__)


def open_file(filename, mode="r"):
    if filename.endswith(".gz"):
        if mode == "w":
            return gzip.open(filename, "wt")
        return gzip.open(filename, "rt")
    else:
        if mode == "w":
            return click.open_file(filename, "w")
        return click.open_file(filename, "r")
