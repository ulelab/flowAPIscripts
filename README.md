# flowAPIscripts
flowrunanalysis.py

A CLI utility to launch CLIP pipeline runs on Flow.bio
 projects.
It authenticates to the Flow API, retrieves project samples, applies simple filters, binds reference files, and submits executions.

Requirements

Python 3.9+

Packages: requests, numpy

A Flow.bio account with project + pipeline access

Installation

Clone this repository and install requirements:

git clone https://github.com/ulelab/flowAPIscripts.git
cd flowAPIscripts
pip install -r requirements.txt   # if you maintain a requirements file


Make the script executable:

chmod +x flowrunanalysis.py

Usage

Basic command structure:

python3 flowrunanalysis.py --pid PROJECT_ID [options]

Required arguments

--pid / --PID : Flow.bio project ID (string)

Optional arguments

--filter sample_name "<glob>" : Filter project samples by name using shell-style wildcards (*, ?).
Example: --filter sample_name "*A"

-n / --num-chunks : Split selected samples into N execution batches using numpy.array_split.
Example: -n 11 runs 110 samples as 11 executions of ~10 each.

--dry-run : Resolve all inputs and print planned submissions without actually running them.

--verbose : Enable DEBUG-level logging.

--log-level : Explicit log level (DEBUG, INFO, WARNING, ERROR).

Examples

Run CLIP pipeline on all samples in project:

python3 flowrunanalysis.py --pid 602442133516131481


Run only samples with names ending in “A”, in 11 batches, dry-run only:

python3 flowrunanalysis.py --pid 602442133516131481 \
  --filter sample_name "*A" \
  -n 11 \
  --dry-run

Features

Authentication: Prompts for Flow.bio username/password, retrieves token.

Pagination: Fetches all project samples via /projects/{pid}/samples until exhausted.

Filtering: Client-side sample name globbing (fnmatch).

Execution batching: Splits samples into N groups for multiple submissions.

Reference binding: Ensures all required genome/index files from prep execution are attached.

Logging: Configurable via flags; execution URLs printed at the end.

Notes

Currently hard-wired to the CLIP pipeline (IDs + version baked in).

Group defaults to sample name; replicate defaults to 1.

Interactive confirmation is requested before the first submission unless --dry-run is used.