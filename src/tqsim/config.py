# This code is part of TQSim.
#
# (C) Copyright Constantine Quantum Technologies, 2022.
#
# This code is licensed under the Apache License, Version 2.0. You may
# obtain a copy of this license in the LICENSE.txt file in the root directory
# of this source tree or at http://www.apache.org/licenses/LICENSE-2.0.
#
# Any modifications or derivative works of this code must retain this
# copyright notice, and modified files need to carry a notice indicating
# that they have been altered from the originals.

import os

PROGRAM_NAME = "tqsim"
CONFIG_PATH = os.path.join(os.path.expanduser("~"), f".{PROGRAM_NAME}")
STORE_PATH = os.path.join(CONFIG_PATH, "store")


def clean_config_path() -> None:
    """Cleans up the configuration directory by removing all files and subdirectories."""
    if os.path.exists(CONFIG_PATH):
        for root, dirs, files in os.walk(CONFIG_PATH, topdown=False):
            for name in files:
                os.remove(os.path.join(root, name))
            for name in dirs:
                os.rmdir(os.path.join(root, name))
        os.rmdir(CONFIG_PATH)
