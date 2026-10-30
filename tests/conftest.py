import os
import sys
from pathlib import Path

import pytest

# Get the absolute path of the project root directory
project_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# Add the project root directory to Python path
sys.path.insert(0, project_root)


@pytest.fixture
def station_model(tmp_path):
    from npzd import NPZDModel

    config = (Path(project_root) / "config.toml").read_text().split("[FORCING]")[0]
    config = config.replace("flag_stn = 1", "flag_stn = 2")
    config_path = tmp_path / "station.toml"
    config_path.write_text(config + '[FORCING]\nsource = "station"\n')
    return NPZDModel(config_path)
