import os

# Default port number to listen on. Can be overriden by passing a port number
# as the first command line argument, e.g. `python tornadoserver.py 1234`
PORT = 7777

# Directory to store previously computed transfer functions, spectra etc. in
DATABASE_DIR = "cache"

# Maximum number of thread pool workers (only required for multi-user usage)
MAX_THREADPOOL_WORKERS = 8

# Path of colormap directory relative to the static directory from which
# tornado serves static files
COLORMAP_PATH = os.path.join("images", "colormaps")

# number of sample points for the transfer function that is displayed
# in the client
TRANSFER_FUNCTION_CLIENT_SAMPLES = 400

# number of sample points for the matter spectrum that is displayed
# in the client per decade
MATTER_SPECTRUM_CLIENT_SAMPLES_PER_DECADE = 40
