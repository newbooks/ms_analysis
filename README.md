# Microstate Analysis Library
## ms_analysis

MCCE step 4 can record all the micorstates Monte Carlo sampling has accepted. Analyzing these accepted states reveals how conformational changes of sites are affected by each other.
This tool offers functions to extract microstates information from the recorded file.

## Mechanism
This library has been rewritten to minimize memory use. The microstates are accessed via index numbes and the real microstate information is saved in an off-line database.

## Data structure and functions
Microstate class:
    Microstate.state:   microstate, an encoded string to represent a list of selected conformers that define the microstate
    Microstate.E:       microstate energy
    Microstate.crg:     microstate net charge (of selected conformers)
    Microstate.counts:  counts of acceptance in MC

MC class:  # A monte Carlo run
    MC.T:               Temperature
    MC.pH:              pH
    MC.Eh:              Eh
    MC.fixedconf:       a list of fixed conformers
    MC.microstates:     A dictionary that uses index numbers as keys and Microstate data structure as value.
        {   index: Microstate, 
            index: Microstate,
            ...
        }

Conformer class:    # conformer info as recorded in head3.lst
    Conformer.iconf:    index, starting from 0
    Conformer.confid:   conformer ID
    Conformer.resid:    residue ID
    Conformer.crg:      charge
    Conformer.load():   read a head3.lst line to load conformer info




## Examples of how to use
from ms_analysis import *

