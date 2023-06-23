# -*- coding: utf-8 -*-
  
# Libraries
import argparse
import gzip
import os

import pandas as pd
import numpy as np

## Given tokenized sequenced and locations that have been identified as either a phage or bacteria, match back to the original sequence.
