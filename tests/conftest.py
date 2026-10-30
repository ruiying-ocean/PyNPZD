import os
import sys
from pathlib import Path

import pytest

# Get the absolute path of the project root directory
project_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# Add the project root directory to Python path
sys.path.insert(0, project_root)


@pytest.fixture
def station_model(monkeypatch):
    from pynpzd import NPZDModel

    monkeypatch.chdir(project_root)
    return NPZDModel(Path(__file__).parent / "fixtures/empower/config.toml")
