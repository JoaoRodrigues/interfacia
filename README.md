## Description
Web server utility to visualize the interfacial energetics of molecular complexes, e.g. protein-protein interactions. **The web server is still in active development and is not suitable for any sort of research.**

## Requirements
The energetics of the interface are calculated using CNS (http://cns-online.org/v1.3) and the OPLS-X force field as implemented in HADDOCK (http://haddocking.org). Please request an executable or the source code of CNS and install it locally to run the web server.

Python library requirements are listed in the `requirements.txt` file, which can be fed to `pip` for an easy installation.

## Example
An example use case (analysis of the E2A-HPR crystal structure) is available under /results/example.

## Usage
````bash
git clone https://github.com/JoaoRodrigues/interfacia.git
cd interfacia
python run.py
````

The web server will be accessible at http://localhost:5000 and an example result page at http://localhost:5000/results/example
