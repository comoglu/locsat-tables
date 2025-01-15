#!/usr/bin/env python3

"""
Seismic Travel Time Table Generator
---------------------------------
Generates travel time tables for various seismic phases using different velocity models.
Supports IASP91 and AK135 models through SeisComP's libtau implementation.

Enhanced features:
- Multiple modes for local, regional, and teleseismic events
- Configurable through command line arguments
- Improved error handling and logging
- Better code organization and maintainability
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

# Use SC to compute the travel times
try:
    from seiscomp.seismology import TravelTimeTable
except ImportError:
    logging.error("SeisComP's Python library not found. Please ensure seiscomp is properly installed.")
    sys.exit(1)

@dataclass
class Config:
    """Configuration class to store all parameters."""
    model: str = "iasp91"  # Velocity model (iasp91 or ak135)
    prefix: str = "iasp91"  # Prefix for output filenames
    conrad_depth: float = 20.0  # Conrad discontinuity depth
    moho_depth: float = 35.0  # Mohorovičić discontinuity depth
    combine_phases: bool = True  # Combine related phases into single table
    extended_format: bool = False  # Use extended output format
    extra_comments: bool = False  # Add extra debugging comments
    output_dir: str = "tables"  # Directory for output files
    mode: str = "default"  # Mode: default, local, regional, or custom
    max_depth: float = 800.0  # Maximum depth in km
    depth_step: float = 5.0  # Depth step in km for shallow depths
    deep_depth_step: float = 50.0  # Depth step in km for deeper depths
    deep_depth_start: float = 100.0  # Depth at which to switch to larger step
    max_distance: float = 180.0  # Maximum distance in degrees
    distance_step: float = 1.0  # Distance step in degrees
    regional_distance_step: float = 0.2  # Distance step for regional distances

    def validate(self) -> None:
        """Validate configuration parameters."""
        if self.model not in ["iasp91", "ak135"]:
            raise ValueError(f"Unsupported model: {self.model}")
        if self.conrad_depth <= 0 or self.moho_depth <= 0:
            raise ValueError("Depths must be positive")
        if self.conrad_depth >= self.moho_depth:
            raise ValueError("Conrad depth must be less than Moho depth")
        if self.max_depth <= 0 or self.depth_step <= 0:
            raise ValueError("Depth values must be positive")
        if self.max_distance <= 0 or self.distance_step <= 0:
            raise ValueError("Distance values must be positive")
        if self.mode not in ["default", "local", "regional", "custom"]:
            raise ValueError("Invalid mode")

@dataclass
class Arrival:
    """Class to store arrival information."""
    phase: str = ""
    time: float = 0.0
    dtdd: float = 0.0  # dt/dΔ: derivative with respect to distance
    dtdh: float = 0.0  # dt/dh: derivative with respect to depth
    
    @classmethod
    def from_traveltime(cls, tt) -> 'Arrival':
        """Create an Arrival from a SeisComP TravelTime object."""
        return cls(
            phase=tt.phase,
            time=tt.time,
            dtdd=tt.dtdd if hasattr(tt, 'dtdd') else 0.0,
            dtdh=tt.dtdh if hasattr(tt, 'dtdh') else 0.0
        )
    
    def __str__(self) -> str:
        return f"{self.phase:<12} {self.time:8.3f} {self.dtdd:8.5f} {self.dtdh:8.5f}"

class TravelTimeCalculator:
    """Handles travel time calculations and table generation."""
    
    def __init__(self, config: Config):
        self.config = config
        self.ttt = TravelTimeTable()
        self.ttt.setModel(config.model)
        self._setup_phase_definitions()
    
    def _setup_phase_definitions(self) -> None:
        """Initialize phase-specific distance and depth ranges."""
        # Setup base distances for different modes
        if self.config.mode == "local":
            self.regional_distances = list(arange(0, 5.1, 0.1))  # Up to 5 degrees, 0.1 degree steps
            max_distance = 5.0
            distance_step = 0.1
        elif self.config.mode == "regional":
            self.regional_distances = list(arange(0, 20.1, 0.2))  # Up to 20 degrees, 0.2 degree steps
            max_distance = 20.0
            distance_step = 0.2
        else:  # default or custom
            self.regional_distances = [
                0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6,
                1.8, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5
            ]
            max_distance = self.config.max_distance
            distance_step = self.config.distance_step

        # Define distance ranges for each phase
        self.distances = {
            "P": self.regional_distances + list(arange(5.0, max_distance + distance_step, distance_step)),
            "Pg": self.regional_distances + list(arange(0, min(10.2, max_distance), self.config.regional_distance_step)),
            "Pb": self.regional_distances + list(arange(0, min(9.5, max_distance), 0.5)),
            "Pn": self.regional_distances + list(arange(0, min(20.5, max_distance), 0.5)),
            "S": self.regional_distances + list(arange(5.0, min(116.0, max_distance), distance_step)),
            "Sb": self.regional_distances + list(arange(0, min(9.5, max_distance), 0.5)),
            "Sg": self.regional_distances + list(arange(0, min(10.2, max_distance), self.config.regional_distance_step)),
            "Sn": self.regional_distances + list(arange(0, min(20.5, max_distance), 0.5)),
        }

        if self.config.mode not in ["local", "regional"]:
            # Add teleseismic phases for default and custom modes
            self.distances.update({
                "pP": arange(20, min(105, max_distance), distance_step),
                "sP": arange(20, min(105, max_distance), distance_step),
                "sS": arange(20, min(105, max_distance), distance_step),
                "PP": arange(30, max_distance + distance_step, distance_step),
                "SS": arange(30, max_distance + distance_step, distance_step),
                "PcP": arange(25, min(61, max_distance), distance_step),
                "ScS": arange(25, min(61, max_distance), distance_step),
                "ScP": arange(25, min(61, max_distance), distance_step),
                "PKP": arange(90, max_distance + distance_step, distance_step),
                "PKPdf": arange(90, max_distance + distance_step, distance_step),
                "PKPab": arange(140, max_distance + distance_step, distance_step),
                "PKPbc": arange(140, min(161, max_distance), distance_step),
                "SKPdf": arange(90, max_distance + distance_step, distance_step),
                "pPKPdf": arange(90, max_distance + distance_step, distance_step),
                "pPKPab": arange(140, max_distance + distance_step, distance_step),
                "pPKPbc": arange(140, min(161, max_distance), distance_step),
                "sPKPdf": arange(90, max_distance + distance_step, distance_step),
                "sPKPab": arange(140, max_distance + distance_step, distance_step),
                "sPKPbc": arange(140, min(161, max_distance), distance_step),
            })

        # Setup depth ranges
        if self.config.mode == "local":
            self.crustal_depths = list(arange(0, 35 + self.config.depth_step, self.config.depth_step))
            max_depth = 35
        elif self.config.mode == "regional":
            self.crustal_depths = list(arange(0, 60 + self.config.depth_step, self.config.depth_step))
            max_depth = 100
        else:
            self.crustal_depths = list(arange(0, 35 + self.config.depth_step, self.config.depth_step))
            max_depth = self.config.max_depth

        # Define depth ranges for different depth regions
        shallow_depths = list(arange(0, min(self.config.deep_depth_start, max_depth), 
                                  self.config.depth_step))
        deep_depths = list(arange(self.config.deep_depth_start, max_depth + self.config.deep_depth_step, 
                                self.config.deep_depth_step))
        self.teleseismic_depths = list(dict.fromkeys(shallow_depths + deep_depths))  # Remove duplicates

        # Set up depths dictionary
        self.depths = {
            "P": self.teleseismic_depths,
            "pP": self.teleseismic_depths,
            "sP": self.teleseismic_depths,
            "Pg": self.crustal_depths,
            "Pb": self.crustal_depths,
            "Pn": self.crustal_depths + [d for d in self.teleseismic_depths if d <= 400],
            "S": self.teleseismic_depths,
            "sS": self.teleseismic_depths,
            "Sb": self.crustal_depths,
            "Sg": self.crustal_depths,
            "Sn": self.crustal_depths + [d for d in self.teleseismic_depths if d <= 400],
        }

        if self.config.mode not in ["local", "regional"]:
            # Add teleseismic phases for default and custom modes
            for phase in ["PP", "SS", "PcP", "ScS", "ScP", "PKP", "PKPdf", "PKPab", 
                        "PKPbc", "SKPdf", "pPKPdf", "pPKPab", "pPKPbc", "sPKPdf", 
                        "sPKPab", "sPKPbc"]:
                self.depths[phase] = self.teleseismic_depths

        # Define distance limits
        self.distance_limits = {
            "Pdiff": (0, min(116, max_distance)),
            "PKPdf": (118, max_distance),
            "Sdiff": (0, min(116, max_distance)),
            "SKSdf": (118, max_distance),
        }

        # Phase combinations remain the same as they define physical relationships
        self.phase_combinations = {
            "P": ["Pg", "Pb", "Pn", "P", "Pdiff", "PKPdf"],
            "PKP": ["PKPab", "PKPbc", "PKPdf"],
            "pP": ["pP", "pPn", "pPdiff"],
            "sP": ["sP", "sPn", "sPdiff"],
            "PP": ["PP", "PnPn"],
            "S": ["S", "Sn", "Sg", "Sb", "Sdiff"],
            "sS": ["sS", "sSg", "sSb", "sSn", "sSdiff"],
            "SS": ["SS", "SnSn"],
            "Pg": ["Pg", "PgPg"],
            "Sg": ["Sg", "SgSg"],
        }

    def compute_travel_times(self, delta: float, depth: float) -> List[Arrival]:
        """
        Compute travel times for given distance and depth.
        
        Args:
            delta: Distance in degrees
            depth: Source depth in km
            
        Returns:
            List of Arrival objects
        """
        try:
            print(f"DEBUG: Computing travel times with delta={delta:.3f}, depth={depth:.3f}")
            
            # Convert inputs to float to ensure correct type
            delta, depth = float(delta), float(depth)
            
            # Use the 6-parameter version as specified in the C++ interface
            print(f"DEBUG: About to call compute(0.0, 0.0, {depth}, 0.0, {delta}, 0.0)")
            arrivals = self.ttt.compute(0.0, 0.0, depth, 0.0, delta, 0.0)
            
            if arrivals:
                print(f"DEBUG: Found {len(arrivals)} arrivals")
                for arr in arrivals:
                    try:
                        print(f"DEBUG: Arrival phase={arr.phase} time={arr.time:.3f} dtdd={arr.dtdd:.5f} dtdh={arr.dtdh:.5f}")
                    except:
                        print(f"DEBUG: Raw arrival: {arr}")
            else:
                print("DEBUG: No arrivals found")
                
            return arrivals
            
        except Exception as e:
            logging.error(f"Error computing travel times for delta={delta}, depth={depth}: {e}")
            return []

    def approximate_pg_time(self, distance: float, depth: float) -> Arrival:
        """Approximate Pg travel time using a simplified velocity model."""
        v = 5.8 + 0.0006 * depth  # Velocity in km/s with slight depth dependence
        time = ((111.195 * distance)**2 + depth**2)**0.5 / v
        return Arrival(phase="Pg", time=time, dtdd=-1, dtdh=-1)

    def approximate_sg_time(self, distance: float, depth: float) -> Arrival:
        """Approximate Sg travel time using a simplified velocity model."""
        v = 3.36 + 0.00037 * depth  # Velocity in km/s with slight depth dependence
        time = ((111.195 * distance)**2 + depth**2)**0.5 / v
        return Arrival(phase="Sg", time=time, dtdd=-1, dtdh=-1)

    def _find_arrival(self, phase: str, distance: float, depth: float) -> Optional[Arrival]:
        """Find the appropriate arrival for given phase, distance and depth."""
        # Special handling for Pg and Sg phases
        if phase == "Pg":
            return self.approximate_pg_time(distance, depth)
        if phase == "Sg":
            return self.approximate_sg_time(distance, depth)

        # Handle depth phases at zero depth
        orig_phase = phase
        if depth == 0 and phase[0] in ["p", "s"]:
            phase = phase[1:]

        arrivals = self.compute_travel_times(distance, depth)
        
        # Check distance limits and phase combinations
        for arr in arrivals:
            if arr.phase in self.distance_limits:
                dmin, dmax = self.distance_limits[arr.phase]
                if not dmin <= distance <= dmax:
                    continue

            if orig_phase in self.phase_combinations:
                if arr.phase in self.phase_combinations[orig_phase]:
                    return Arrival.from_traveltime(arr)
            elif arr.phase == orig_phase:
                return Arrival.from_traveltime(arr)

        return None

    def create_table(self, phase: str) -> str:
        """Create a travel time table for a specific phase."""
        if phase not in self.depths or phase not in self.distances:
            raise ValueError(f"Unknown phase: {phase}")

        lines = []
        
        # Phase name without comments
        lines.append(phase)
        
        # Write depth samples
        depths = self.depths[phase]
        lines.append(f"{len(depths)}    # number of depth samples")
        lines.extend(self._format_number_list(depths))
        
        # Write distance samples
        distances = self.distances[phase]
        lines.append(f"{len(distances)}    # number of distance samples")
        lines.extend(self._format_number_list(distances))

        # Process each depth
        for depth in depths:
            # Only output required depth header
            lines.append(f"# z = {depth} km")
            
            # Process each distance for this depth
            for distance in distances:
                arrival = self._find_arrival(phase, distance, depth)
                if arrival:
                    if depth == 0.0 and distance == 0.0:
                        arrival.time = 0.0
                    # Only output the time value
                    lines.append(f"{arrival.time:12.3f}")
                else:
                    lines.append(f"{-1:12.3f}")

        # LocSAT format doesn't want trailing newline
        return "\n".join(lines)

    @staticmethod
    def _format_number_list(numbers: List[float], width: int = 8) -> List[str]:
        """Format a list of numbers into columns."""
        lines = []
        current_line = []
        
        for num in numbers:
            current_line.append(f"{num:{width}.2f}")
            if len(current_line) == 10:  # 10 numbers per line
                lines.append("".join(current_line))
                current_line = []
        
        if current_line:
            # Pad last line to match width
            while len(current_line) < 10:
                current_line.append(" " * width)
            lines.append("".join(current_line))
            
        return lines

def main():
    """Main function to handle command line interface and program execution."""
    parser = argparse.ArgumentParser(
        description="Generate seismic travel time tables.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    parser.add_argument(
        "--model", 
        choices=["iasp91", "ak135"], 
        default="iasp91",
        help="Velocity model to use"
    )
    
    parser.add_argument(
        "--output-dir", 
        default="tables",
        help="Output directory for travel time tables"
    )
    
    parser.add_argument(
        "--mode",
        choices=["default", "local", "regional", "custom"],
        default="default",
        help="""Mode of operation:
            default: Standard tables for all distances and depths
            local: Optimized for local earthquakes (0-5 degrees, 0-35 km)
            regional: Optimized for regional earthquakes (0-20 degrees, 0-100 km)
            custom: Use custom distance and depth ranges"""
    )
    
    # Custom depth parameters
    parser.add_argument(
        "--max-depth",
        type=float,
        default=800.0,
        help="Maximum depth in km (for custom mode)"
    )
    
    parser.add_argument(
        "--depth-step",
        type=float,
        default=5.0,
        help="Depth step in km for shallow depths (for custom mode)"
    )
    
    parser.add_argument(
        "--deep-depth-step",
        type=float,
        default=50.0,
        help="Depth step in km for deeper depths (for custom mode)"
    )
    
    parser.add_argument(
        "--deep-depth-start",
        type=float,
        default=100.0,
        help="Depth at which to switch to larger depth step (for custom mode)"
    )
    
    # Custom distance parameters
    parser.add_argument(
        "--max-distance",
        type=float,
        default=180.0,
        help="Maximum distance in degrees (for custom mode)"
    )
    
    parser.add_argument(
        "--distance-step",
        type=float,
        default=1.0,
        help="Distance step in degrees for teleseismic distances (for custom mode)"
    )
    
    parser.add_argument(
        "--regional-distance-step",
        type=float,
        default=0.2,
        help="Distance step in degrees for regional distances (for custom mode)"
    )
    
    # Other parameters
    parser.add_argument(
        "--conrad-depth",
        type=float,
        default=20.0,
        help="Depth of Conrad discontinuity in km"
    )
    
    parser.add_argument(
        "--moho-depth",
        type=float,
        default=35.0,
        help="Depth of Mohorovičić discontinuity in km"
    )
    
    parser.add_argument(
        "--combine-phases",
        action="store_true",
        default=True,
        help="Combine related phases into single table"
    )
    
    parser.add_argument(
        "--verbosity",
        choices=['DEBUG', 'INFO', 'WARNING', 'ERROR'],
        default='INFO',
        help="Set the logging level"
    )
    
    parser.add_argument(
        "phases",
        nargs="*",
        help="Specific phases to compute (default: all supported phases)"
    )
    
    args = parser.parse_args()
    
    # Setup logging
    logging.basicConfig(
        level=getattr(logging, args.verbosity),
        format='%(asctime)s - %(levelname)s - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )
    
    # Define default phases based on mode
    if args.mode == "local":
        default_phases = "P Pg Pb Pn S Sg Sb Sn".split()
    elif args.mode == "regional":
        default_phases = "P Pg Pb Pn S Sg Sb Sn pP sP PP SS".split()
    else:
        default_phases = """
            P Pg Pb Pn S Sg Sb Sn
            PP SS
            pP sP sS
            PKP PKPab PKPbc PKPdf pPKPab pPKPbc pPKPdf sPKPab sPKPbc sPKPdf SKPdf
            PcP ScS ScP
        """.strip().split()
    
    phases = args.phases if args.phases else default_phases
    
    try:
        # Initialize configuration
        config = Config(
            model=args.model,
            prefix=args.model,
            conrad_depth=args.conrad_depth,
            moho_depth=args.moho_depth,
            combine_phases=args.combine_phases,
            extended_format=False,  # Always use standard format for LocSAT
            extra_comments=False,   # No extra comments for LocSAT
            output_dir=args.output_dir,
            mode=args.mode,
            max_depth=args.max_depth,
            depth_step=args.depth_step,
            deep_depth_step=args.deep_depth_step,
            deep_depth_start=args.deep_depth_start,
            max_distance=args.max_distance,
            distance_step=args.distance_step,
            regional_distance_step=args.regional_distance_step
        )
        
        # Validate configuration
        config.validate()
        
        # Create output directory
        output_path = Path(config.output_dir)
        output_path.mkdir(parents=True, exist_ok=True)
        
        # Initialize calculator
        calculator = TravelTimeCalculator(config)
        
        total_phases = len(phases)
        logging.info(f"Starting travel time calculations for {total_phases} phases")
        logging.info(f"Using velocity model: {config.model}")
        logging.info(f"Mode: {config.mode}")
        logging.info(f"Output directory: {output_path}")
        
        # Process each phase
        for i, phase in enumerate(phases, 1):
            try:
                logging.info(f"Processing phase {i}/{total_phases}: {phase}")
                table = calculator.create_table(phase)
                output_file = output_path / f"{config.prefix}.{phase}"
                # Write with explicit Unix newlines and no trailing newline
                with open(output_file, 'w', newline='\n') as f:
                    f.write(table)
                logging.info(f"Successfully created table: {output_file}")
                
            except Exception as e:
                logging.error(f"Error processing phase {phase}: {e}")
                if args.verbosity == 'DEBUG':
                    logging.exception("Detailed error information:")
                continue
        
        logging.info("Travel time table generation completed")
        
    except Exception as e:
        logging.error(f"Fatal error: {e}")
        if args.verbosity == 'DEBUG':
            logging.exception("Detailed error information:")
        sys.exit(1)

if __name__ == "__main__":
    main()