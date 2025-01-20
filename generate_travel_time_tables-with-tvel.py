#!/usr/bin/env python3

"""
TVEL-based Travel Time Table Generator
------------------------------------
Generates travel time tables using velocity model files in TVEL format.
Combines functionality from local and global travel time table generators.
"""

import sys
import os
import logging
import argparse
from typing import Dict, List, Tuple, Optional
from dataclasses import dataclass
from pathlib import Path
import numpy as np

@dataclass
class LocalVelocityModel:
    """Class to store velocity model data from TVEL file."""
    depths: np.ndarray
    vp: np.ndarray
    vs: np.ndarray
    density: np.ndarray
    max_depth: float

    @classmethod
    def from_tvel_file(cls, filename: str) -> 'LocalVelocityModel':
        """Create a LocalVelocityModel from a TVEL file."""
        depths, vp, vs, density = [], [], [], []
        
        with open(filename, 'r') as f:
            # Skip header line
            header = f.readline()
            model_name = header.strip()
            logging.info(f"Loading velocity model: {model_name}")
            
            for line in f:
                line = line.strip()
                if not line:  # Skip empty lines
                    continue
                
                try:
                    # Handle fixed-width format with potential extra spaces
                    parts = line.split()
                    if len(parts) != 4:
                        continue
                        
                    depth = float(parts[0])
                    vp_val = float(parts[1])
                    vs_val = float(parts[2])
                    den_val = float(parts[3])
                    
                    # Skip if all values are zero (sometimes used as layer markers)
                    if vp_val == 0 and vs_val == 0 and den_val == 0:
                        continue
                        
                    depths.append(depth)
                    vp.append(vp_val)
                    vs.append(vs_val)
                    density.append(den_val)
                    
                except (ValueError, IndexError) as e:
                    logging.warning(f"Skipping invalid line: {line.strip()}")
        
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
        """Get velocity at a specific depth using linear interpolation."""
        if depth > self.max_depth:
            return None
            
        idx = np.searchsorted(self.depths, depth)
        if idx == 0:
            velocity = self.vp[0] if wave_type.upper() == 'P' else self.vs[0]
            return velocity if velocity != 0 else None
        elif idx == len(self.depths):
            velocity = self.vp[-1] if wave_type.upper() == 'P' else self.vs[-1]
            return velocity if velocity != 0 else None
            
        # Linear interpolation
        d0, d1 = self.depths[idx-1], self.depths[idx]
        if wave_type.upper() == 'P':
            v0, v1 = self.vp[idx-1], self.vp[idx]
        else:
            v0, v1 = self.vs[idx-1], self.vs[idx]
            
        return v0 + (v1 - v0) * (depth - d0) / (d1 - d0)

class TravelTimeCalculator:
    """Handles travel time calculations using velocity model."""
    
    def __init__(self, velocity_model: LocalVelocityModel):
        self.velocity_model = velocity_model
        
    def calculate_time(self, phase: str, depth: float, distance: float) -> Optional[float]:
        """Calculate travel time for given phase, depth and distance."""
        if depth > self.velocity_model.max_depth:
            return None
            
        # Convert distance to kilometers
        distance_km = distance * 111.195
        
        # Determine wave type from phase
        if phase.upper().startswith('P'):
            velocity = self.velocity_model.get_velocity(depth, 'P')
        elif phase.upper().startswith('S'):
            velocity = self.velocity_model.get_velocity(depth, 'S')
        else:
            return None
            
        if velocity is None:
            return None
            
        # Simple straight-line approximation
        hypotenuse = np.sqrt(distance_km**2 + depth**2)
        return hypotenuse / velocity

def format_number_list(numbers: List[float], width: int = 8) -> List[str]:
    """Format a list of numbers into columns."""
    lines = []
    current_line = []
    
    for num in numbers:
        current_line.append(f"{num:{width}.2f}")
        if len(current_line) == 10:  # 10 numbers per line
            lines.append("".join(current_line))
            current_line = []
    
    if current_line:
        while len(current_line) < 10:
            current_line.append(" " * width)
        lines.append("".join(current_line))
        
    return lines

def create_travel_time_table(
    velocity_model: LocalVelocityModel,
    phase: str,
    depth_samples: List[float] = [0, 5, 15, 30, 40, 50, 75, 100, 150, 200, 300, 400, 500, 600, 800],
    distance_range: Tuple[float, float, float] = (0, 180, 1)
) -> str:
    """Create a travel time table for a specific phase in the required format."""
    calculator = TravelTimeCalculator(velocity_model)
    lines = []
    
    # Write phase name with comment
    lines.append(f"{phase} # travel-time (and amplitude) tables")
    
    # Write depth samples
    lines.append(f"{len(depth_samples)}    # number of depth samples")
    # Format depth samples in rows of 10
    for i in range(0, len(depth_samples), 10):
        line_samples = depth_samples[i:i+10]
        line = "".join(f"{depth:8.2f}" for depth in line_samples)
        lines.append(line)
    
    # Generate distance samples
    start, end, step = distance_range
    distances = [i * step for i in range(int((end - start) / step) + 1)]
    
    # Write distance samples
    lines.append(f"{len(distances)}    # number of distance samples")
    # Format distance samples in rows of 10
    for i in range(0, len(distances), 10):
        line_samples = distances[i:i+10]
        line = "".join(f"{dist:8.2f}" for dist in line_samples)
        lines.append(line)
    
    # Calculate travel times for each depth
    for depth in depth_samples:
        lines.append(f"# Travel-time/amplitude for z = {depth:8.2f}")
        for distance in distances:
            time = calculator.calculate_time(phase, depth, distance)
            if time is not None:
                lines.append(f"{time:12.3f}")
            else:
                lines.append(f"{-1:12.3f}")
    
    return "\n".join(lines)

def main():
    parser = argparse.ArgumentParser(description='Generate travel time tables from TVEL files')
    
    # Required arguments
    parser.add_argument('tvel_file', help='Input TVEL velocity model file')
    
    # Optional arguments
    parser.add_argument('--output-dir', default='tables',
                      help='Output directory for travel time tables')
    
    parser.add_argument('--phases', type=str, default='P,S,Pg,Sg',
                      help='Comma-separated list of phases (default: P,S,Pg,Sg)')
    
    parser.add_argument('--depth-samples', type=str, 
                      default='0,5,15,30,40,50,75,100,150,200,300,400,500,600,800',
                      help='Comma-separated list of depth samples in km')
    
    parser.add_argument('--distance-range', type=str, default='0,180,1',
                      help='Distance range as start,end,step in degrees (default: 0,180,1)')
    
    parser.add_argument('--verbosity', choices=['DEBUG', 'INFO', 'WARNING', 'ERROR'],
                      default='INFO', help='Set logging level')
    
    args = parser.parse_args()
    
    args = parser.parse_args()
    
    try:
        # Setup logging
        logging.basicConfig(
            level=getattr(logging, args.verbosity),
            format='%(asctime)s - %(levelname)s - %(message)s'
        )

        # Parse parameters
        depth_samples = [float(x) for x in args.depth_samples.split(',')]
        dstart, dend, dstep = map(float, args.distance_range.split(','))
        phases = [p.strip() for p in args.phases.split(',')]

        # Create output directory
        output_path = Path(args.output_dir)
        output_path.mkdir(parents=True, exist_ok=True)
        
        # Load velocity model
        model = LocalVelocityModel.from_tvel_file(args.tvel_file)
        model_name = Path(args.tvel_file).stem
        
        # Generate tables for each phase
        for phase in phases:
            output_file = output_path / f"{model_name}.{phase}"
            
            logging.info(f"Generating table for phase {phase}")
            table = create_travel_time_table(
                model,
                phase,
                depth_samples=depth_samples,
                distance_range=(dstart, dend, dstep)
            )
            
            with open(output_file, 'w', newline='\n') as f:
                f.write(table)
            
            print(f"Generated {output_file}")
            
    except Exception as e:
        logging.error(f"Error generating travel time tables: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()