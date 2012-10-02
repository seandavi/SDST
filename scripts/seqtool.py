#!/usr/bin/env python
from seqtools.main import parser
from seqtools.commands import *

args = parser.parse_args()
args.func(args)
