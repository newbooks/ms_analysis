#!/usr/bin/env python
# Charge microstate analysis tool sets
# This library offers the following functions
# * Read microstate file and convert to a list of charge microstates
# * Group microstates by residue charge state
# * Group/rank microstates by population
# * Bin microstates
# * Find microstates within an energy band
# * Distance score of two microstate groups

import math
import numpy as np
import sys

ph2Kcal = 1.364
Kcal2kT = 1.688


class Microstate:
    def __init__(self, state, E, count):
        self.state = ""   # selected conformers of a microstate
        self.E = 0
        self.crg = 0
        self.count = 0


class Conformer:
    def __init__(self):
        self.iconf = 0
        self.confid = ""
        self.resid = ""
        self.crg = 0.0

    def load(self, line):
        fields = line.split()
        self.iconf = int(fields[0]) - 1
        self.confid = fields[1]
        self.resid = self.confid[:3]+self.confid[5:11]
        self.crg = float(fields[4])

class Residue:
    def __init__(self):
        self.resid = ""
        self.conformers = []

class MC:
    def __init__(self):
        self.T = 298.15
        self.pH = 7.00
        self.Eh = 0.00
        self.method = ""
        self.counts = 0
        self.E = 0.0
        self.conformers = []
        self.fixedconfs = []                # fixed conformers
        self.free_residues = []             # a list of conformer groups that make up free residues
        self.ires_by_iconf = {}             # index of free residue by index of conf
        self.microstates = []
        lines = open("head3.lst").readlines()
        for line in lines[1:]:
            if len(line) > 80:
                conf = Conformer()
                conf.load(line)
                self.conformers.append(conf)

    def readms(self, fname):
        found_nres = False
        found_mc = False
        newmc = False

        f = open(fname)

        # read the header part
        # mc condition
        fields = []
        n_lines = 0
        while len(fields) != 3 and n_lines < 10:
            line = f.readline()
            line = line.split("#")[0]  #remove comment
            fields = line.split(",")
            n_lines += 1
        if n_lines >= 10:
            print("Expect MC condition line like \"T:298.15,pH:5.00,eH:0.00\" in the first 10 lines")
            sys.exit()
        for field in fields:
            key, value = field.split(":")
            key = key.strip()
            value = float(value.strip())
            if key.upper() == "T":
                self.T = value
            elif key.upper() == "PH":
                self.pH = value
            elif key.upper() == "EH":
                self.Eh = value

        # method
        fields = []
        while len(fields) != 2:
            line = f.readline()
            line = line.split("#")[0]  #remove comment
            fields = line.split(":")
        if fields[0].strip() == "METHOD":
            self.method = fields[1].strip()
        else:
            print("Expect line of METHOD record")
            sys.exit()

        # fixed conformer, skip
        fields = []
        while len(fields) != 2:
            line = f.readline()
            line = line.split("#")[0]  #remove comment
            fields = line.split(":")

        # free residues


        # read MC microstates


        f.close()


def readheadlst(fname):
    conformers = []
    lines = open(fname).readlines()

    for line in lines[1:]:
        if len(line) > 80:
            conf = Conformer()
            conf.load(line)
            conformers.append(conf)

    return conformers


if __name__ == "__main__":
    msfile = "ms_out/pH5eH0ms.txt"

    mc = MC()
    mc.readms("ms_out/pH4eH0ms.txt")
