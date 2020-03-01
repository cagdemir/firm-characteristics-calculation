
import pandas as pd
import numpy as np

# 0. write necessary functions - ranking function
# 1. get data paths
# 2. prepare global frame and train frame
# 3. read return, pass it to global frame, shift it by -1, prepare return mask, flatten it and attach to train frame
# 4. read var, pass it to global frame, apply return mask, rank it, multiply w/ macroeconomcis data and attach them to train frame
