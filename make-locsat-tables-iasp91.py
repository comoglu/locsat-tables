#!/usr/bin/env python

import sys
import os
from typing import List, Dict, Tuple, Union
from dataclasses import dataclass
import numpy as np
from seiscomp.seismology import TravelTimeTable
import argparse
import logging
from tqdm import tqdm

# Configuration
MODEL = "iasp91"  # Can be either "iasp91" or "ak135"
PREFIX = MODEL  # for filename generation
CONRAD_DEPTH = 20.0
MOHO_DEPTH = 35.0
COMBINE_PHASES = True
EXTENDED_FORMAT = False
EXTRA_COMMENTS = False

# Setup logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

@dataclass
class Arrival:
    phase: str
    time: float
    dtdd: float
    dtdh: float

    def __str__(self):
        return f"{self.phase:<12} {self.time:8.3f} {self.dtdd:8.5f} {self.dtdh:8.5f}"

class TravelTimeCalculator:
    def __init__(self, model: str):
        self.ttt = TravelTimeTable()
        self.ttt.setModel(model)

    def compute(self, delta: float, depth: float) -> List[Arrival]:
        arrivals = self.ttt.compute(0, 0, depth, 0, delta, 0, 0)
        return [Arrival(arr.phase, arr.time, arr.dtdd, arr.dtdh) for arr in arrivals]

def print_as_10_columns(data: Union[List[float], np.ndarray]) -> str:
    data = np.array(data)
    return '\n'.join(' '.join(f'{x:8.2f}' for x in row) for row in data.reshape(-1, 10))

def find_arrival(arrivals: List[Arrival], phase: str, distance: float) -> Union[Arrival, None]:
    for arr in arrivals:
        if arr.phase in distanceLimits:
            dmin, dmax = distanceLimits[arr.phase]
            if not dmin <= distance <= dmax:
                continue

        if phase in phaseCombinations:
            if arr.phase in phaseCombinations[phase]:
                return arr
        elif arr.phase == phase:
            return arr
    return None

def approximate_pg(depth: float, distance: float) -> Arrival:
    v = 5.8 + 0.0006 * depth
    time = ((111.195 * distance) ** 2 + depth ** 2) ** 0.5 / v
    return Arrival("Pg", time, -1, -1)

def approximate_sg(depth: float, distance: float) -> Arrival:
    v = 3.36 + 0.00037 * depth
    time = ((111.195 * distance) ** 2 + depth ** 2) ** 0.5 / v
    return Arrival("Sg", time, -1, -1)

def create_table(phase: str, calculator: TravelTimeCalculator) -> str:
    output = f"{'#' if not EXTENDED_FORMAT else ''} {phase} travel-time tables\n"
    output += f"{len(depths[phase])}    # number of depth samples\n"
    output += print_as_10_columns(depths[phase])
    output += f"{len(distances[phase])}    # number of distance samples\n"
    output += print_as_10_columns(distances[phase])

    for z in tqdm(depths[phase], desc=f"Processing {phase}"):
        ph = phase[1:] if z == 0 and phase[0] in ["p", "s"] else phase

        if EXTENDED_FORMAT:
            output += f"# z = {z} km\n"
        else:
            output += f"# Travel time for z = {z} km\n"

        for d in distances[phase]:
            arrivals = calculator.compute(d, z)
            arr = find_arrival(arrivals, ph, d)

            if phase == "Pg":
                arr = approximate_pg(z, d)
            elif phase == "Sg":
                arr = approximate_sg(z, d)

            if EXTRA_COMMENTS:
                output += f"# depth: {z:8.3f}    distance: {d:.3f}\n"

            if arr:
                if z == 0.0 and d == 0.0:
                    arr.time = 0.0
                if EXTENDED_FORMAT:
                    output += f"{arr.time:12.3f} {arr.dtdd:12.6f} {arr.dtdh:12.6f} {arr.phase}\n"
                else:
                    output += f"{arr.time:12.3f}\n"
            else:
                if EXTENDED_FORMAT:
                    output += f"{-1:12.3f} {-1:12.6f} {-1:12.6f}\n"
                else:
                    output += f"{-1:12.3f}\n"

    return output

# Define regional_distances, distances, crustal_depths, teleseismic_depths, depths, distanceLimits, and phaseCombinations here
# (These are large dictionaries and lists, so I'm omitting them for brevity)

def parse_arguments():
    parser = argparse.ArgumentParser(description="Generate seismic travel-time tables.")
    parser.add_argument("--model", choices=["iasp91", "ak135"], default=MODEL, help="Velocity model to use")
    parser.add_argument("--combine-phases", action="store_true", default=COMBINE_PHASES, help="Combine phases like Pg, Pn, P, Pdiff, PKPdf into a single table")
    parser.add_argument("--extended-format", action="store_true", default=EXTENDED_FORMAT, help="Generate extended format output")
    parser.add_argument("--extra-comments", action="store_true", default=EXTRA_COMMENTS, help="Add extra comments for debugging")
    parser.add_argument("phases", nargs="*", help="Phases to generate tables for")
    return parser.parse_args()

def main():
    args = parse_arguments()
    global MODEL, COMBINE_PHASES, EXTENDED_FORMAT, EXTRA_COMMENTS
    MODEL = args.model
    COMBINE_PHASES = args.combine_phases
    EXTENDED_FORMAT = args.extended_format
    EXTRA_COMMENTS = args.extra_comments

    calculator = TravelTimeCalculator(MODEL)
    phase_set = args.phases if args.phases else defaultPhaseSet.strip().split()

    try:
        os.makedirs("tables", exist_ok=True)
    except OSError as e:
        logger.error(f"Error creating 'tables' directory: {e}")
        sys.exit(1)

    for phase in phase_set:
        try:
            table = create_table(phase, calculator)
            filename = os.path.join("tables", f"{PREFIX}.{phase}")
            with open(filename, "w") as f:
                f.write(table)
            logger.info(f"Generated table for {phase}")
        except Exception as e:
            logger.error(f"Error generating table for {phase}: {e}")

if __name__ == "__main__":
    main()