## Description
Web server utility to visualize the interfacial energetics of molecular complexes, e.g. protein-protein interactions. The energetics of the interface are calculated using CNS (http://cns-online.org/v1.3) and the OPLS-X force field as implemented in HADDOCK (http://haddocking.org). **The web server is still in active development and is not suitable for any sort of research.**

## Requirements
### Software
* [Crystallography and NMR System (CNS) v1.3](http://cns-online.org/v1.3)
* [Redis](http://redis.io/)
* Python 2.7+

## Installation Instructions
After obtaining CNS v1.3 and redis, use `pip` to install the Python library requirements listed in the `requirements.txt` file.

```bash
git clone https://github.com/JoaoRodrigues/interfacia.git
cd interfacia
pip install -r requirements.txt
```

## Running the web server
```bash
# Start the redis server
redis-server --daemonize yes
# Start the Celery queue
celery worker -f celery.log -D --quiet -A app.celery_app 
# Start the Flask web server
python run.py

# To stop the web server:
# 1. Kill the running Flask process
# 2. Kill the Celery queue process
ps auxww | grep 'celery worker' | awk '{print $2}' | xargs kill
# 3. Stop the redis-server daemon
redis-cli shutdown
```

The web server will be accessible at http://localhost:5000 and an example result page at http://localhost:5000/results/example

## Example
An example use case (analysis of the Barnase-Barstar crystal structure) is available under /results/example.
