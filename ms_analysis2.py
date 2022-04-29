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
import tracemalloc

ph2Kcal = 1.364
Kcal2kT = 1.688


class Microstate:
    def __init__(self, state, E, count):
        self.stateid = " ".join([str(x) for x in state])
        self.E = E
        self.count = count

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
        self.microstates = []               # a list of microstates
        self.microstates_by_id = {}
        lines = open("head3.lst").readlines()
        for line in lines[1:]:
            if len(line) > 80:
                conf = Conformer()
                conf.load(line)
                self.conformers.append(conf)
        self.allms = []

    def readms(self, fname):
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
        fields = []
        while len(fields) != 2:
            line = f.readline()
            line = line.split("#")[0]  #remove comment
            fields = line.split(":")

        n_res = int(fields[0])
        self.free_residues = [[int(xx) for xx in x.strip().split()] for x in fields[1].strip(" ;\n").split(";")]

        if len(self.free_residues) != n_res:
            print("The number of free residues don't match.")
            print(line)
            sys.exit()

        for ires in range(len(self.free_residues)):
            res = self.free_residues[ires]
            for iconf in res:
                self.ires_by_iconf[iconf] = ires

        # read MC microstates
        newmc = False
        found_mc = False
        self.microstates_by_id.clear()
        while True:
            line = f.readline()
            if not line:
                break

            if line.find("MC:") == 0:   # found a MC start
                newmc = True
                found_mc = True
                continue
            elif newmc:
                f1, f2 = line.split(":")
                current_state = [int(c) for c in f2.split()]
                newmc = False
                continue
            elif found_mc:
                fields = line.split(",")
                if len(fields) >= 3:
                    state_e = float(fields[0])
                    count = int(fields[1])
                    flipped = [int(c) for c in fields[2].split()]

                    for ic in flipped:
                        ir = self.ires_by_iconf[ic]
                        current_state[ir] = ic

                    #print(flipped, current_state)
                    ms = Microstate(current_state, state_e, count)
                    if ms.stateid in self.microstates_by_id:
                        self.microstates_by_id[ms.stateid].count += ms.count
                    else:
                        self.microstates_by_id[ms.stateid] = ms

        f.close()

        # convert microstates to a list
        self.microstates = [item[1] for item in self.microstates_by_id.items()]
        self.allms = range(len(self.microstates))


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
    tracemalloc.start()
    mc = MC()
    mc.readms(msfile)
    print("Loaded ms", tracemalloc.get_traced_memory())

    ms_odd = []
    ms_even = []
    for i in range(len(mc.microstates)):
        if i % 2:
            ms_even.append(mc.microstates[i])
        else:
            ms_odd.append(mc.microstates[i])

    print("divide ms", tracemalloc.get_traced_memory())

    tracemalloc.stop()
    print(len(ms_odd), len(ms_even))