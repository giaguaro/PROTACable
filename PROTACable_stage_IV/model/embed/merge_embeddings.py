"""
------------------------------------------------------------------------------------------
Copyright (C) 2023 Hazem Mslati. All Rights Reserved.

This script is provided "AS IS" without any warranty.
Unauthorized copying, modification, or distribution is prohibited.
------------------------------------------------------------------------------------------
"""

import glob
import numpy as np
import joblib
import os


featurized_files_esm = glob.glob('*featurized_esm.pkl')
new_files = glob.glob('*lframe.pkl')

featurized_files_esm = featurized_files_esm[0]
new_files = new_files[0]
        
featurized_file_esm = joblib.load(featurized_files_esm)
new_file = joblib.load(new_files)

new_file.append(featurized_file_esm[-2:])
new_file.append(featurized_file_esm[-1:])


# Save combined data
joblib.dump(new_file, new_files, compress=3)
