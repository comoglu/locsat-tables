# Seismic Travel-time Table Generator

## Overview

This Python script generates travel-time tables for various seismic phases using either the IASP91 or AK135 velocity models. It's an improved version of the original script by Joachim Saul, featuring enhanced error handling, logging, command-line arguments, and overall code structure.

## Features

- Generate travel-time tables for multiple seismic phases
- Support for IASP91 and AK135 velocity models
- Command-line interface for easy configuration
- Progress bar for long-running computations
- Comprehensive error handling and logging
- Modular code structure for easy maintenance and extension

## Requirements

- Python 3.7+
- NumPy
- tqdm
- seiscomp (make sure this is installed and properly configured)

## Installation

1. Clone this repository:
   ```
   git clone https://github.com/your-username/seismic-traveltime-generator.git
   cd seismic-traveltime-generator
   ```

2. Install the required packages:
   ```
   pip install numpy tqdm
   ```

   Note: The `seiscomp` package needs to be installed separately according to its own installation instructions.

## Usage

Run the script from the command line with various options:

```
python scttgen.py [OPTIONS] [PHASES...]
```

### Options:

- `--model {iasp91,ak135}`: Choose the velocity model (default: iasp91)
- `--combine-phases`: Combine phases like Pg, Pn, P, Pdiff, PKPdf into a single table
- `--extended-format`: Generate extended format output
- `--extra-comments`: Add extra comments for debugging

### Examples:

1. Generate tables for P and S phases using the IASP91 model:
   ```
   python scttgen.py P S
   ```

2. Use the AK135 model with combined phases and extended format:
   ```
   python scttgen.py --model ak135 --combine-phases --extended-format P S PKP
   ```

3. Generate all default phases with extra debugging comments:
   ```
   python scttgen.py --extra-comments
   ```

## Output

The script generates travel-time tables in the `tables/` directory. Each file is named according to the pattern `{model}.{phase}`, e.g., `iasp91.P` for the P-phase table using the IASP91 model.

## Contributing

Contributions to improve the script are welcome. Please feel free to submit issues or pull requests.

## License

This project is licensed under the [insert appropriate license, e.g., MIT License, GPL, etc.]. See the LICENSE file for details.

## Acknowledgements

This script is based on the original work by Joachim Saul. The improvements were made to enhance usability, maintainability, and extensibility.

