# RUN) python ena_file_rename.py [directory/path]

import os
import sys

# current dir
dir_path = sys.argv[1]
file_name = os.listdir(dir_path)

# file rename
for name in file_name:
    # raw_name
    src = os.path.join(dir_path, name)

    split_name = name.split("_")
    dst = split_name[0] + "_R" + split_name[1]
    dst = os.path.join(dir_path, dst)

    os.renames(src, dst)

    









