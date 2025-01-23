# LocSAT Travel Time Table Generator

A collection of Python scripts for generating travel time tables used by LocSAT (Location based on Satellite data) for seismic event location. This toolkit supports various velocity models including IASP91, AK135, and custom local velocity models.

## Scripts Overview

- `generate_travel_time_tables.py`: Main script for generating tables using IASP91/AK135 models
- `generate_local_tables.py`: Generates tables using local velocity models
- `generate_travel_time_tables-with-tvel.py`: Creates tables from TVEL format velocity models

## Requirements

- Python 3.6+
- SeisComP (with Python bindings)
- NumPy

## Installation

1. Clone the repository:
```bash
git clone https://github.com/comoglu/locsat-tables.git
cd locsat-tables
```

2. Ensure SeisComP is installed with Python bindings.

## Usage

### Global Models (IASP91/AK135)

```bash
python generate_travel_time_tables.py --model iasp91 --output-dir tables [OPTIONS]
```

Options:
- `--model`: Choose between "iasp91" or "ak135" (default: iasp91)
- `--mode`: Operating mode ["default", "local", "regional", "custom"]
- `--output-dir`: Output directory for tables
- `--max-depth`: Maximum depth in km (default: 800.0)
- `--depth-step`: Depth step for shallow depths (default: 5.0)
- `--distance-step`: Distance step in degrees (default: 1.0)

### Local Velocity Models

```bash
python generate_local_tables.py velocity_model.txt --output-dir tables [OPTIONS]
```

Options:
- `--reference`: Reference model for depths beyond local model
- `--prefix`: Prefix for output files
- `--output-dir`: Output directory for tables

### TVEL Format Models

```bash
python generate_travel_time_tables-with-tvel.py model.tvel --output-dir tables [OPTIONS]
```

Options:
- `--phases`: Comma-separated list of phases (default: P,S,Pg,Sg)
- `--depth-samples`: Comma-separated depth values in km
- `--distance-range`: Distance range as start,end,step in degrees

## Velocity Model Formats

### Local Velocity Model Format
```
depth(km)  vp(km/s)  vs(km/s)  density(g/cm³)
0.0        5.8       3.36      2.7
...
```

### TVEL Format
```
ModelName
depth(km)  vp(km/s)  vs(km/s)  density(g/cm³)
0.0        5.8       3.36      2.7
...
```

## Supported Phases

- Local/Regional: P, Pg, Pb, Pn, S, Sg, Sb, Sn
- Teleseismic: P, S, PP, SS, PKP, SKS, etc.

## Output Format

The generated tables follow the LocSAT format:
```
phase_name
number_of_depth_samples
depth_samples
number_of_distance_samples
distance_samples
travel_times
```

## Error Handling

- The scripts include comprehensive error handling and logging
- Use `--verbosity DEBUG` for detailed error information
- Invalid or missing data points are marked with -1 in the output

## Contributing

1. Fork the repository
2. Create a feature branch
3. Submit a pull request

## License

MIT License