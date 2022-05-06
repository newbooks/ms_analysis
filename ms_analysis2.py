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
import zlib

ph2Kcal = 1.364
Kcal2kT = 1.688


class Microstate:
    def __init__(self, state, E, count):
        self.stateid = zlib.compress(" ".join([str(x) for x in state]).encode())
        self.E = E
        self.count = count

    def state(self):
        return [int(i) for i in zlib.decompress(self.stateid).decode().split()]

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

class Free_Res:
    def __init__(self):
        self.resid = ""
        self.charges = []     # a list of charge choice

class Charge_Microstate:
    def __init__(self, crg_stateid, total_E, count):
        self.crg_stateid = crg_stateid
        self.average_E = 0
        self.total_E = total_E
        self.count = count

    def state(self):
        return [int(i) for i in zlib.decompress(self.crg_stateid).decode().split()]

def readheadlst(fname):
    conformers = []
    lines = open(fname).readlines()

    for line in lines[1:]:
        if len(line) > 80:
            conf = Conformer()
            conf.load(line)
            conformers.append(conf)

    return conformers


class MC:
    def __init__(self):
        self.T = 298.15
        self.pH = 7.00
        self.Eh = 0.00
        self.method = ""
        self.counts = 0
        self.E = 0.0
        self.conformers = []
        self.iconf_by_confname = {}
        self.fixedconfs = []                # fixed conformers
        self.free_residues = []             # a list of conformer groups that make up free residues
        self.free_residue_names = []        # free residue names
        self.ires_by_iconf = {}             # index of free residue by index of conf
        self.microstates = []               # a list of microstates
        self.microstates_by_id = {}
        lines = open("head3.lst").readlines()
        iconf = 0
        for line in lines[1:]:
            if len(line) > 80:
                conf = Conformer()
                conf.load(line)
                self.conformers.append(conf)
                self.iconf_by_confname[conf.confid] = iconf
                iconf += 1

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
        self.free_residue_names = [x[0] for x in self.free_residues]

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
                    self.counts += ms.count

        f.close()

        # convert microstates to a list
        self.microstates = [item[1] for item in self.microstates_by_id.items()]

    def get_occ(self, microstates):
        conf_occ = [0 for _ in range(len(self.conformers))]
        total_counts = 0
        for ms in microstates:
            total_counts += ms.count
            for iconf in ms.state():
                conf_occ[iconf] += ms.count
        for ic in range(len(conf_occ)):
            conf_occ[ic] = conf_occ[ic]/total_counts

        return conf_occ

    def select_by_conformer(self, microstates, conformer_in=[]):
        "Select microstate if confomer is in the list AND energy is in the range. Return all if the list is empty."
        selected = []
        unselected = []
        if conformer_in:
            iconf_in = set([self.iconf_by_confname[confid] for confid in conformer_in])
        else:
            return [], microstates

        for ms in microstates:
            state = set(ms.state())
            if state & iconf_in:
                selected.append(ms)
            else:
                unselected.append(ms)

        return selected, unselected

    def select_by_energy(self, microstates, energy_in=[]):
        "Select microstate if energy is in the list AND energy is in the range. Return all if the list is empty."
        selected = []
        unselected = []
        if energy_in:
            energy_in.sort()
        else:
            return [], microstates

        for ms in microstates:
            if energy_in[0] <= ms.E < energy_in[1]:
                selected.append(ms)
            else:
                unselected.append(ms)

        return selected, unselected

    def convert_to_charge_ms(self):
        charge_microstates = []
        charge_ms_by_id = {}   #
        for ms in self.microstates:
            current_crg_state = [round(self.conformers[ic].crg) for ic in ms.state()]
            crg_stateid = zlib.compress(" ".join([str(x) for x in current_crg_state]).encode())
            crg_ms = Charge_Microstate(crg_stateid, ms.E*ms.count, ms.count)
            if crg_stateid in charge_ms_by_id:
                charge_ms_by_id[crg_stateid].count += crg_ms.count
                charge_ms_by_id[crg_stateid].total_E += crg_ms.total_E
            else:
                charge_ms_by_id[crg_stateid] = crg_ms
        for crg_stateid in charge_ms_by_id.keys():
            crg_ms = charge_ms_by_id[crg_stateid]
            crg_ms.average_E = crg_ms.total_E / crg_ms.count
            charge_microstates.append(crg_ms)
        del(charge_ms_by_id)
        return charge_microstates

def get_erange(microstates):
    "return energy range of the microstates"
    emin = microstates[0].E
    emax = microstates[0].E
    for ms in microstates[1:]:
        if emin > ms.E:
            emin = ms.E
        if emax < ms.E:
            emax = ms.E
    return emin, emax


def bin_mscounts_total(microstates, nbins=100, erange=[]):
    "Return two lists, one as energy range and one as counts of total counts."
    if erange:   # use the range if range arrary is provided
        erange.sort()
        counts = [0 for _ in erange]
        for ms in microstates:
            ibin = -1
            for ie in range(len(erange)-1, -1, -1):
                if ms.E > erange[ie]:
                    ibin = ie
                    break
            if ibin >= 0:
                counts[ibin] += ms.count
    else:
        lowest_E, highest_E = get_erange(microstates)

        E_range = highest_E - lowest_E + 0.01
        bin_size = E_range / nbins
        counts = [0 for _ in range(nbins)]
        for ms in microstates:
            ibin = int((ms.E - lowest_E) / bin_size)
            counts[ibin] += ms.count
        erange = [lowest_E + i*bin_size for i in range(nbins)]

    return erange, counts

def bin_mscounts_unique(microstates, nbins=100, erange=[]):
    "Return two lists, one as energy range and one as counts of to counts."

    if erange:   # use the range if range arrary is provided
        erange.sort()
        counts = [0 for _ in erange]
        for ms in microstates:
            ibin = -1
            for ie in range(len(erange)-1, -1, -1):
                if ms.E > erange[ie]:
                    ibin = ie
                    break
            if ibin >= 0:
                counts[ibin] += 1
    else:
        lowest_E, highest_E = get_erange(microstates)

        E_range = highest_E - lowest_E + 0.01
        bin_size = E_range / nbins
        counts = [0 for _ in range(nbins)]
        for ms in microstates:
            ibin = int((ms.E - lowest_E) / bin_size)
            counts[ibin] += 1
        erange = [lowest_E + i*bin_size for i in range(nbins)]

    return erange, counts


def get_count(microstates):
    "Calculate the microstate count"
    count = 0
    for ms in microstates:
        count += ms.count
    return count


def average_e(microstates):
    t = 0.0
    c = 0
    for ms in microstates:
        t += ms.E * ms.counts
        c += ms.counts
    return t/c

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

if __name__ == "__main__":
    msfile = "ms_out/pH5eH0ms.txt"
    #tracemalloc.start()
    mc = MC()
    mc.readms(msfile)
    #print("Loaded ms", tracemalloc.get_traced_memory())


    # # Example 1: Bin microstates based on energy
    # erange, total_counts = bin_mscounts_total(mc.microstates)
    # erange, uniq_counts = bin_mscounts_unique(mc.microstates)
    # for i in range(len(erange)):
    #     print("%8.3f %6d %6d" % (erange[i], total_counts[i], uniq_counts[i]))
    # #print("bin ms", tracemalloc.get_traced_memory())
    #
    #
    # # Example 2: When GLU35 is ionized, what residues change conformation?
    # glu_charged_confs = ["GLU-1A0035_011", "GLU-1A0035_012", "GLU-1A0035_013", "GLU-1A0035_011"]
    # glu_charged_ms, glu_neutral_ms = mc.select_by_conformer(mc.microstates, conformer_in=glu_charged_confs)
    # conf_occ_glu_charged = mc.get_occ(glu_charged_ms)
    # conf_occ_glu_neutral = mc.get_occ(glu_neutral_ms)
    # for res in mc.free_residues:
    #     resid = mc.conformers[res[0]].resid
    #     prob1 = [conf_occ_glu_neutral[ic] for ic in res]
    #     prob2 = [conf_occ_glu_charged[ic] for ic in res]
    #     d = bhata_distance(prob1, prob2)
    #     print("%s, d= %.3f" % (resid, d))
    #     for ic in res:
    #         print("%s %6.3f %6.3f" % (mc.conformers[ic].confid, conf_occ_glu_neutral[ic], conf_occ_glu_charged[ic]))
    #     print()
    # #print("Grouping ms", tracemalloc.get_traced_memory())

    # Example 3: Which charge microstate is the most dominant?
    charge_microstates = mc.convert_to_charge_ms()
    charge_microstates.sort(key=lambda x: x.count)
    count = 0
    for crg_ms in charge_microstates:
        count += crg_ms.count
        print(crg_ms.state(), crg_ms.count, crg_ms.average_E)
    print("%d charge microstates" % len(charge_microstates))
    print("%d total microstates" % count)
    #print("charge microstates", tracemalloc.get_traced_memory())


    #tracemalloc.stop()
