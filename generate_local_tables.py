#!/usr/bin/env python3

"""
Local Velocity Model Travel Time Table Generator
---------------------------------------------
Generates travel time tables for various seismic phases using local velocity models,
with the option to merge with global models (IASP91/AK135) for deeper structures.
"""

import sys
import os
import logging
import argparse
from typing import Dict, List, Tuple, Optional, Union
from dataclasses import dataclass
from pathlib import Path
import numpy as np
from numpy import arange, array, concatenate

try:
    from seiscomp.seismology import TravelTimeTable
except ImportError:
    logging.error("SeisComP's Python library not found. Please ensure seiscomp is properly installed.")
    sys.exit(1)

@dataclass
class LocalVelocityModel:
    """Class to store local velocity model data."""
    depths: np.ndarray
    vp: np.ndarray
    vs: np.ndarray
    density: np.ndarray
    max_depth: float

    @classmethod
    def from_file(cls, filename: str) -> 'LocalVelocityModel':
        """Create a LocalVelocityModel from a file."""
        depths, vp, vs, density = [], [], [], []
        
        with open(filename, 'r') as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith('#'):
                    continue
                    
                values = [float(x) for x in line.split()]
                if len(values) == 4:
                    depths.append(values[0])
                    vp.append(values[1])
                    vs.append(values[2])
                    density.append(values[3])
        
        depths = np.array(depths)
        max_depth = np.max(depths)
        
        return cls(
            depths=np.array(depths),
            vp=np.array(vp),
            vs=np.array(vs),
            density=np.array(density),
            max_depth=max_depth
        )

    def get_velocity(self, depth: float, wave_type: str) -> float:
        """Get velocity at a specific depth."""
        if depth > self.max_depth:
            return None
            
        layer_index = np.searchsorted(self.depths, depth) - 1
        if layer_index < 0:
            layer_index = 0
            
        if wave_type.upper() == 'P':
            return self.vp[layer_index]
        elif wave_type.upper() == 'S':
            return self.vs[layer_index]
        return None

class LocalTravelTimeCalculator:
    """Handles travel time calculations using local velocity model."""
    
    def __init__(self, local_model: LocalVelocityModel, reference_model: str = "ak135"):
        self.local_model = local_model
        self.reference_model = reference_model
        self.ttt = TravelTimeTable()
        self.ttt.setModel(reference_model)
        
    def calculate_local_time(self, phase: str, depth: float, distance: float) -> Optional[float]:
        """Calculate travel time using local velocity model."""
        min_dist, max_dist, max_depth, velocity_type = self._get_phase_characteristics(phase)
        
        if distance < min_dist or distance > max_dist:
            return None
            
        distance_km = distance * 111.19
        
        if depth > self.local_model.max_depth:
            # Use reference model for depths beyond local model
            try:
                arrivals = self.ttt.compute(0.0, 0.0, depth, 0.0, distance, 0.0)
                for arr in arrivals:
                    if arr.phase == phase:
                        return arr.time
                return None
            except Exception as e:
                logging.error(f"Error computing travel time with reference model: {e}")
                return None
        
        # Get velocity based on wave type
        if velocity_type.startswith('P'):
            velocity = self.local_model.get_velocity(depth, 'P')
        elif velocity_type.startswith('S'):
            velocity = self.local_model.get_velocity(depth, 'S')
        else:
            velocity = None
            
        if velocity is None:
            return None
            
        # Simple straight-line approximation
        hypotenuse = np.sqrt(distance_km**2 + depth**2)
        return hypotenuse / velocity
    
    @staticmethod
    def _get_phase_characteristics(phase: str) -> Tuple[float, float, float, str]:
        """Get the characteristics for a specific phase."""
        characteristics = {
            'P': (0, 180, 800, 'P'),
            'S': (0, 180, 800, 'S'),
            'Pg': (0, 10, 30, 'P'),
            'Pb': (0, 10, 30, 'P'),
            'Sg': (0, 10, 30, 'S'),
            'Sb': (0, 10, 30, 'S'),
            'Pn': (2, 15, 50, 'P'),
            'Sn': (2, 15, 50, 'S'),
        }
        return characteristics.get(phase, (0, 180, 800, 'P'))

def main():
    parser = argparse.ArgumentParser(description='Generate travel time tables using local velocity model')
    parser.add_argument('velocity_model', help='Local velocity model file')
    parser.add_argument('--reference', choices=['ak135', 'iasp91'], default='ak135',
                      help='Reference model for depths beyond local model')
    parser.add_argument('--output-dir', default='tables',
                      help='Output directory for travel time tables')
    parser.add_argument('--prefix', help='Prefix for output files')
    
    args = parser.parse_args()
    
    try:
        # Create output directory
        output_path = Path(args.output_dir)
        output_path.mkdir(parents=True, exist_ok=True)
        
        # Read local velocity model
        local_model = LocalVelocityModel.from_file(args.velocity_model)
        
        # Set prefix based on velocity model file if not provided
        prefix = args.prefix or Path(args.velocity_model).stem
        
        # Initialize calculator
        calculator = LocalTravelTimeCalculator(local_model, args.reference)
        
        # Define phases for local/regional events
        phases = ['P', 'S', 'Pg', 'Pb', 'Pn', 'Sg', 'Sb', 'Sn']
        
        for phase in phases:
            output_file = output_path / f"{prefix}_{args.reference}.{phase}"
            
            # Generate travel time table
            depth_samples = np.arange(0, min(local_model.max_depth + 5, 800), 5)
            distance_samples = np.arange(0, 20.1, 0.2)  # Up to 20 degrees for regional
            
            table_lines = []
            table_lines.append(phase)
            table_lines.append(f"{len(depth_samples)}    # number of depth samples")
            
            # Format depth samples
            for i in range(0, len(depth_samples), 10):
                line = "    " + "    ".join(f"{d:7.2f}" for d in depth_samples[i:i+10])
                table_lines.append(line)
                
            table_lines.append(f"{len(distance_samples)}    # number of distance samples")
            
            # Format distance samples
            for i in range(0, len(distance_samples), 10):
                line = "    " + "    ".join(f"{d:7.2f}" for d in distance_samples[i:i+10])
                table_lines.append(line)
                
            # Calculate travel times
            for depth in depth_samples:
                table_lines.append(f"# z = {depth:8.2f} km")
                for distance in distance_samples:
                    time = calculator.calculate_local_time(phase, depth, distance)
                    if time is not None:
                        table_lines.append(f"{time:12.3f}")
                    else:
                        table_lines.append(f"{-1:12.3f}")
            
            # Write table to file
            with open(output_file, 'w', newline='\n') as f:
                f.write('\n'.join(table_lines))
            
            print(f"Generated {output_file}")
            
    except Exception as e:
        logging.error(f"Error generating travel time tables: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()