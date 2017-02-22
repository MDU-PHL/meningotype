#!/usr/bin/env python
# Global functions

# Use modern print function from python 3.x
from __future__ import print_function

# Modules
import os
import sys

# Log a message to stderr
def msg(*args, **kwargs):
	print(*args, file=sys.stderr, **kwargs)

# Log an error to stderr and quit with non-zero error code
def err(*args, **kwargs):
	msg(*args, **kwargs)
	sys.exit(1);
