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
import matplotlib.pyplot as plt
import numpy as np
import sys

ph2Kcal = 1.364
Kcal2kT = 1.688


class Microstate:
    def __init__(self, state, E, count):
        self.state = state   # a list of charges on free residues
        self.E_sum = 0
        self.E = 0           # Average energy
        self.count = count


class Conformer:
    def __init__(self):
        self.iconf = 0
        self.confid = ""
        self.resid = ""
        self.occ = 0.0
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
        self.residues = []
        self.ionizable_residues = []      # a list of conformer groups that make up ionizable residues
        self.ires_by_confname = {}
        self.ires_by_iconf = {}
        self.ires_by_resname = {}
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
        fields = []
        while len(fields) != 2:
            line = f.readline()
            line = line.split("#")[0]  #remove comment
            fields = line.split(":")
        res_fields = fields[1].strip(" \n;").split(";")
        for res_field in res_fields:
            iconfs = [int(i) for i in res_field.split()]
            charges = []
            for iconf in iconfs:
                charges.append(round(self.conformers[iconf].crg))
            if len(set(charges)) > 1:     # more than one charge, ionizable
                self.ionizable_residues.append(list(set(charges)))
                for iconf in iconfs:
                    self.ires_by_iconf[iconf] = len(self.ionizable_residues) - 1
            print(self.ionizable_residues)

        # read MC microstates


        f.close()

        #
        #
        #     line = f.readline()
        #     if line.find("#N_FREE residues:") == 0:
        #         found_nres = True
        #         continue
        #     elif found_nres:
        #         fields = line.split(":")
        #         fields = fields[1].split(";")
        #         for f in fields:
        #             if f.strip():
        #                 confs = f.split()
        #                 free_res.append([int(c) for c in confs])
        #         for i_res in range(len(free_res)):
        #             for iconf in free_res[i_res]:
        #                 iconf2res[iconf] = i_res
        #         found_nres = False
        #     elif line.find("MC:") == 0:   # ms starts
        #         found_mc = True
        #         newmc = True
        #         continue
        #     elif newmc:
        #         f1, f2 = line.split(":")
        #         current_state = [int(c) for c in f2.split()]
        #         newmc = False
        #         continue
        #     elif found_mc:
        #         fields = line.split(",")
        #         if len(fields) >= 3:
        #             state_e = float(fields[0])
        #             count = int(fields[1])
        #             flipped = [int(c) for c in fields[2].split()]
        #
        #             for ic in flipped:
        #                 ir = iconf2res[ic]
        #                 current_state[ir] = ic
        #
        #             ms = Microstate(list(current_state), state_e, count)
        #             microstates.append(ms)
        #
        # return microstates, free_res, iconf2res


def readheadlst(fname):
    conformers = []
    lines = open(fname).readlines()

    for line in lines[1:]:
        if len(line) > 80:
            conf = Conformer()
            conf.load(line)
            conformers.append(conf)

    return conformers


def get_occ(microstates):
    "Return a dict index by iconf, value as counts and occupancy"
    conf_counts = {}
    n_states = 0
    for ms in microstates:
        for iconf in ms.state:
            if iconf in conf_counts:
                conf_counts[iconf] += ms.count
            else:
                conf_counts[iconf] = ms.count
        n_states += ms.count

    conf_occ = {}
    for key in conf_counts.keys():
        conf_occ[key] = (conf_counts[key], conf_counts[key]/n_states)

    return conf_occ


def groupms_byconf(microstates, confs):
    # One of ic in confs has to be in microstate for microstate to be matched
    matched = []
    unmatched = []
    for ms in microstates:
        match = False
        for ic in confs:
            if ic in ms.state:
                match = True
                break
        if match:
            matched.append(ms)
        else:
            unmatched.append(ms)

    return matched, unmatched


def groupms_byenergy(microstates, bounds):
    # One of ic in confs has to be in microstate for microstate to be matched
    matched = []
    unmatched = []
    for ms in microstates:
        low_bound = min(bounds)
        high_bound = max(bounds)

        if low_bound <= ms.E <= high_bound:
            matched.append(ms)
        else:
            unmatched.append(ms)

    return matched, unmatched


def bhata_distance(prob1, prob2):
    if len(prob1) != len(prob2):
        d = 1.0e10  # Max possible value set to this
    else:
        t = 0.0
        for i in range(len(prob1)):
            t += math.sqrt(prob1[i]*prob2[i])
        bc = t
        if bc < math.exp(-100):
            d = 100.0
        else:
            d = -math.log(bc)

    return d






def ms_counts_total(microstates, nbins = 100, erange = []):
    if erange:   # use the range if range arrary is provided
        counts = [0 for _ in erange]
        for ms in microstates:
            ibin = -1
            for ie in range(len(erange)-1):
                if erange[ie] <= ms.E < erange[ie+1]:
                    ibin = ie
                    break
            if ibin >= 0:
                counts[ibin] += ms.count
    else:
        lowest_E = highest_E = microstates[0].E

        for ms in microstates[1:]:
            if lowest_E > ms.E:
                lowest_E = ms.E
            if highest_E < ms.E:
                highest_E = ms.E

        E_range = highest_E - lowest_E + 0.01
        bin_size = E_range / nbins
        counts = [0 for _ in range(nbins)]
        for ms in microstates:
            ibin = int((ms.E - lowest_E) / bin_size)
            counts[ibin] += ms.count
        erange = [lowest_E + i*bin_size for i in range(nbins)]

    return erange, counts


def ms_counts_uniq(microstates, nbins = 100, erange = []):
    # reduce microstates to unique set
    ms_ids = set()

    if erange:   # use the range if range arrary is provided
        counts = [0 for _ in erange]
        for ms in microstates:
            if str(ms.state) in ms_ids:   # identical ms already counted
                continue
            ibin = -1
            for ie in range(len(erange)-1):
                if erange[ie] <= ms.E < erange[ie+1]:
                    ibin = ie
                    break
            if ibin >= 0:
                counts[ibin] += 1
                ms_ids.add(str(ms.state))
            else:
                ms_ids.add(str(ms.state))   # ignore out off range point next time

    else:
        lowest_E = highest_E = microstates[0].E

        for ms in microstates[1:]:
            if lowest_E > ms.E:
                lowest_E = ms.E
            if highest_E < ms.E:
                highest_E = ms.E

        E_range = highest_E - lowest_E + 0.01
        bin_size = E_range / nbins
        counts = [0 for _ in range(nbins)]
        for ms in microstates:
            ibin = int((ms.E - lowest_E) / bin_size)
            if ibin >= 0 and str(ms.state) not in ms_ids:
                counts[ibin] += 1
                ms_ids.add(str(ms.state))
        erange = [lowest_E + i*bin_size for i in range(nbins)]

    return erange, counts

if __name__ == "__main__":
    msfile = "ms_out/pH5eH0ms.txt"

    mc = MC()
    mc.readms(msfile)