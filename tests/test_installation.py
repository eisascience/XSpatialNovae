"""Test model download and cache configuration."""

import os
import sys
from pathlib import Path
from unittest import mock

import pytest


def test_download_models_script_exists():
    """Test that the download_models.py script exists and is executable."""
    script_path = Path(__file__).parent.parent / "scripts" / "download_models.py"
    assert script_path.exists(), f"Script not found at {script_path}"
    
    # Check if it's executable (on Unix-like systems)
    if sys.platform != "win32":
        assert os.access(script_path, os.X_OK), "Script is not executable"


def test_download_models_help():
    """Test that the download_models.py script can show help."""
    import subprocess
    
    script_path = Path(__file__).parent.parent / "scripts" / "download_models.py"
    result = subprocess.run(
        [sys.executable, str(script_path), "--help"],
        capture_output=True,
        text=True
    )
    
    assert result.returncode == 0, f"Script failed: {result.stderr}"
    assert "Download and cache Novae models" in result.stdout
    assert "--novae-model-id" in result.stdout
    assert "--cache-dir" in result.stdout


def test_get_model_cache_dir_default():
    """Test get_model_cache_dir returns default location."""
    # Import the function directly from the module file to avoid full package import
    import importlib.util
    spec = importlib.util.spec_from_file_location(
        "novae_runner", 
        Path(__file__).parent.parent / "novae_seurat_gui" / "modeling" / "novae_runner.py"
    )
    novae_runner = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(novae_runner)
    
    get_model_cache_dir = novae_runner.get_model_cache_dir
    
    # Clear relevant env vars
    env_backup = {}
    for var in ["HF_HOME", "TRANSFORMERS_CACHE", "HUGGINGFACE_HUB_CACHE"]:
        env_backup[var] = os.environ.pop(var, None)
    
    try:
        cache_dir = get_model_cache_dir()
        expected = Path.home() / ".cache" / "huggingface" / "hub"
        assert cache_dir == expected, f"Expected {expected}, got {cache_dir}"
    finally:
        # Restore env vars
        for var, value in env_backup.items():
            if value is not None:
                os.environ[var] = value


def test_get_model_cache_dir_hf_home():
    """Test get_model_cache_dir respects HF_HOME."""
    # Import the function directly from the module file to avoid full package import
    import importlib.util
    import tempfile
    spec = importlib.util.spec_from_file_location(
        "novae_runner", 
        Path(__file__).parent.parent / "novae_seurat_gui" / "modeling" / "novae_runner.py"
    )
    novae_runner = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(novae_runner)
    
    get_model_cache_dir = novae_runner.get_model_cache_dir
    
    test_path = os.path.join(tempfile.gettempdir(), "test_hf_home")
    with mock.patch.dict(os.environ, {"HF_HOME": test_path}, clear=False):
        cache_dir = get_model_cache_dir()
        assert cache_dir == Path(test_path) / "hub"


def test_get_model_cache_dir_transformers_cache():
    """Test get_model_cache_dir respects TRANSFORMERS_CACHE."""
    # Import the function directly from the module file to avoid full package import
    import importlib.util
    import tempfile
    spec = importlib.util.spec_from_file_location(
        "novae_runner", 
        Path(__file__).parent.parent / "novae_seurat_gui" / "modeling" / "novae_runner.py"
    )
    novae_runner = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(novae_runner)
    
    get_model_cache_dir = novae_runner.get_model_cache_dir
    
    test_path = os.path.join(tempfile.gettempdir(), "test_transformers")
    # Clear HF_HOME to test priority
    env = {"TRANSFORMERS_CACHE": test_path}
    if "HF_HOME" in os.environ:
        env["HF_HOME"] = None
    
    with mock.patch.dict(os.environ, env, clear=False):
        # Ensure HF_HOME is not set
        if "HF_HOME" in os.environ:
            del os.environ["HF_HOME"]
        
        cache_dir = get_model_cache_dir()
        assert cache_dir == Path(test_path)


def test_get_model_cache_dir_huggingface_hub_cache():
    """Test get_model_cache_dir respects HUGGINGFACE_HUB_CACHE."""
    # Import the function directly from the module file to avoid full package import
    import importlib.util
    import tempfile
    spec = importlib.util.spec_from_file_location(
        "novae_runner", 
        Path(__file__).parent.parent / "novae_seurat_gui" / "modeling" / "novae_runner.py"
    )
    novae_runner = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(novae_runner)
    
    get_model_cache_dir = novae_runner.get_model_cache_dir
    
    test_path = os.path.join(tempfile.gettempdir(), "test_hf_hub")
    # Clear HF_HOME and TRANSFORMERS_CACHE to test priority
    with mock.patch.dict(os.environ, 
                         {"HUGGINGFACE_HUB_CACHE": test_path}, 
                         clear=False):
        # Ensure higher-priority vars are not set
        if "HF_HOME" in os.environ:
            del os.environ["HF_HOME"]
        if "TRANSFORMERS_CACHE" in os.environ:
            del os.environ["TRANSFORMERS_CACHE"]
        
        cache_dir = get_model_cache_dir()
        assert cache_dir == Path(test_path)


def test_run_novae_zeroshot_has_cache_dir_param():
    """Test that run_novae_zeroshot accepts cache_dir parameter."""
    import inspect
    # Import the function directly from the module file to avoid full package import
    import importlib.util
    spec = importlib.util.spec_from_file_location(
        "novae_runner", 
        Path(__file__).parent.parent / "novae_seurat_gui" / "modeling" / "novae_runner.py"
    )
    novae_runner = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(novae_runner)
    
    run_novae_zeroshot = novae_runner.run_novae_zeroshot
    
    sig = inspect.signature(run_novae_zeroshot)
    assert "cache_dir" in sig.parameters, "cache_dir parameter not found"
    
    # Check it's optional (has default)
    param = sig.parameters["cache_dir"]
    assert param.default is None, "cache_dir should default to None"


def test_requirements_files_exist():
    """Test that all requirements files exist."""
    repo_root = Path(__file__).parent.parent
    
    required_files = [
        "requirements.txt",
        "requirements-uv.txt",
        "requirements-histology.txt",
    ]
    
    for filename in required_files:
        filepath = repo_root / filename
        assert filepath.exists(), f"Missing requirements file: {filename}"
        
        # Check it's not empty
        content = filepath.read_text()
        assert len(content) > 0, f"Requirements file is empty: {filename}"


def test_requirements_include_huggingface_hub():
    """Test that huggingface-hub is in requirements."""
    repo_root = Path(__file__).parent.parent
    
    for req_file in ["requirements.txt", "requirements-uv.txt"]:
        content = (repo_root / req_file).read_text()
        assert "huggingface-hub" in content, f"huggingface-hub not in {req_file}"


def test_pyproject_python_version():
    """Test that pyproject.toml has correct Python version requirements."""
    import toml
    
    repo_root = Path(__file__).parent.parent
    pyproject = toml.load(repo_root / "pyproject.toml")
    
    python_version = pyproject["project"]["requires-python"]
    assert ">=3.11" in python_version, "Should require Python >= 3.11"
    assert "<3.14" in python_version, "Should limit Python < 3.14"


def test_pyproject_has_histology_extras():
    """Test that pyproject.toml defines histology optional dependencies."""
    import toml
    
    repo_root = Path(__file__).parent.parent
    pyproject = toml.load(repo_root / "pyproject.toml")
    
    optional_deps = pyproject["project"]["optional-dependencies"]
    assert "histology" in optional_deps, "histology extras not defined"
    
    # Check histology includes expected packages
    histology_deps = optional_deps["histology"]
    assert any("openslide" in dep for dep in histology_deps), "openslide not in histology extras"


def test_pyproject_includes_huggingface_hub():
    """Test that pyproject.toml includes huggingface-hub in dependencies."""
    import toml
    
    repo_root = Path(__file__).parent.parent
    pyproject = toml.load(repo_root / "pyproject.toml")
    
    deps = pyproject["project"]["dependencies"]
    assert any("huggingface-hub" in dep for dep in deps), "huggingface-hub not in dependencies"
