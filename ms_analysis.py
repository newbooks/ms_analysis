#!/usr/bin/env python
import math
import matplotlib.pyplot as plt
import numpy as np

ph2Kcal = 1.364
Kcal2kT = 1.688


class Microstate:
    def __init__(self, state, E, count):
        self.state = state
        self.E = E
        self.count = count


class Conformer:
    def __init__(self):
        self.iconf = 0
        self.confid = ""
        self.resid = ""
        self.oocc = 0.0
        self.occ = 0.0
        self.crg = 0.0

    def load(self, line):
        fields = line.split()
        self.iconf = int(fields[0]) - 1
        self.confid = fields[1]
        self.resid = self.confid[:3]+self.confid[5:11]
        self.crg = float(fields[4])


def readms(fname):
    current_state = []     # An array of conformer index number as a microstate
    microstates = []       # An array of Microstate objects (conformer index numbers, energy, and count).
    free_res = []          # A list if lists (residues) which contain conformer index numbers
    iconf2res = {}         # A dictionary for looking up from conformer index to residue index

    lines = open(fname).readlines()

    found_nres = False
    found_mc = False
    newmc = False
    for line in lines:
        if line.find("#N_FREE residues:") == 0:
            found_nres = True
            continue
        elif found_nres:
            fields = line.split(":")
            fields = fields[1].split(";")
            for f in fields:
                if f.strip():
                    confs = f.split()
                    free_res.append([int(c) for c in confs])
            for i_res in range(len(free_res)):
                for iconf in free_res[i_res]:
                    iconf2res[iconf] = i_res
            found_nres = False
        elif line.find("MC:") == 0:   # ms starts
            fields = line.split(":")
            which_mc = fields[1].strip()
            found_mc = True
            newmc = True
            continue
        elif newmc:
            f1, f2 = line.split(":")
            current_state = [int(c) for c in f2.split()]
            newmc = False
            continue
        elif found_mc:
            if which_mc in MC_Segments:
                fields = line.split(",")
                if len(fields) >= 3:
                    state_e = float(fields[0])
                    count = int(fields[1])
                    flipped = [int(c) for c in fields[2].split()]

                    for ic in flipped:
                        ir = iconf2res[ic]
                        current_state[ir] = ic

                    ms = Microstate(list(current_state), state_e, count)
                    microstates.append(ms)

    return microstates, free_res, iconf2res


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


def autolabel(rects, ax):
    """Attach a text label above each bar in *rects*, displaying its height."""
    for rect in rects:
        height = rect.get_height()
        ax.annotate('%.3f' % height,
                    xy=(rect.get_x() + rect.get_width() / 2, height),
                    xytext=(0, 3),  # 3 points vertical offset
                    textcoords="offset points",
                    ha='center', va='bottom')


def plot_prob(p1, p2, d):
    x = np.arange(len(p1))
    width = 0.35


    fig, ax = plt.subplots(figsize=(8, 5))
    rects1 = ax.bar(x - width/2, p1, width, label="group1")
    rects2 = ax.bar(x + width/2, p2, width, label="group2")
    ax.set_ylabel("occ")
    ax.set_title("d=%.3f" %d)
    ax.legend()
    autolabel(rects1, ax)
    autolabel(rects2, ax)

    plt.show()
    return


def average_e(microstates):
    t = 0.0
    c = 0
    for ms in microstates:
        t += ms.E * ms.counts
        c += ms.counts
    return t/c


def e2occ(e):
    "Energy (kcal/mol) to occupancy"
    occ_g = 1.0
    occ_x = math.exp(-e*Kcal2kT)
    return occ_x/(occ_g+occ_x)


def occ2e(occ):
    if 1-occ < 0.0001:
        nocc =0.0001

    nocc = occ/(1.0-occ)  # normalized occ, assume ground state G = 0
    return -math.log(nocc)/Kcal2kT


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
