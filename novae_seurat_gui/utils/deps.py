"""Dependency checking and management utilities."""

import importlib
import logging
import shutil
import sys
from typing import Dict, List, Optional, Tuple

logger = logging.getLogger(__name__)


class MissingDependency(Exception):
    """Exception raised when a required dependency is missing."""

    def __init__(self, package_name: str, install_hint: str):
        self.package_name = package_name
        self.install_hint = install_hint
        super().__init__(f"Missing dependency: {package_name}\n{install_hint}")


def is_uv_available() -> bool:
    """
    Check if uv is available in the environment.
    
    Returns
    -------
    bool
        True if uv is available, False otherwise.
    """
    return shutil.which("uv") is not None


def is_in_virtualenv() -> bool:
    """
    Check if running inside a virtual environment.
    
    Returns
    -------
    bool
        True if in a virtual environment, False otherwise.
    """
    # Check for common venv indicators
    return (
        hasattr(sys, 'real_prefix') or  # virtualenv
        (hasattr(sys, 'base_prefix') and sys.base_prefix != sys.prefix) or  # venv
        'VIRTUAL_ENV' in os.environ  # Environment variable check
    )


def get_install_hint(
    package_name: str,
    pip_package: Optional[str] = None,
    prefer_conda: bool = False,
) -> str:
    """
    Generate an install hint for a missing package.
    
    Parameters
    ----------
    package_name : str
        The Python import name of the package.
    pip_package : str, optional
        The pip package name if different from import name.
        If None, uses package_name.
    prefer_conda : bool
        If True, suggest conda first (only for non-venv environments).
        Default is False.
    
    Returns
    -------
    str
        Install hint message.
    """
    pip_pkg = pip_package or package_name
    in_venv = is_in_virtualenv()
    uv_available = is_uv_available()
    
    # Build install command suggestions
    commands = []
    
    # For virtual environments, never suggest --user or conda
    if in_venv:
        if uv_available:
            commands.append(f"uv pip install {pip_pkg}")
        commands.append(f"pip install {pip_pkg}")
    else:
        # Not in venv - suggest conda if preferred and not using uv
        if prefer_conda and not uv_available:
            commands.append(f"conda install {pip_pkg}")
        
        if uv_available:
            commands.append(f"uv pip install {pip_pkg}")
        
        commands.append(f"pip install {pip_pkg}")
    
    if len(commands) == 1:
        return f"Install with: {commands[0]}"
    else:
        cmd_list = "\n  - ".join(commands)
        return f"Install with one of:\n  - {cmd_list}"


def require_package(
    import_name: str,
    pip_package: Optional[str] = None,
    prefer_conda: bool = False,
) -> None:
    """
    Require a package to be installed, raising MissingDependency if not found.
    
    Parameters
    ----------
    import_name : str
        The Python import name of the package (e.g., 'igraph', 'skmisc').
    pip_package : str, optional
        The pip/conda package name if different from import name.
        Examples: 'python-igraph' for 'igraph', 'scikit-misc' for 'skmisc'.
    prefer_conda : bool
        If True and not in a venv, suggest conda installation first.
    
    Raises
    ------
    MissingDependency
        If the package cannot be imported.
    """
    try:
        importlib.import_module(import_name)
        logger.debug(f"Package '{import_name}' is available")
    except ImportError:
        hint = get_install_hint(import_name, pip_package, prefer_conda)
        logger.error(f"Missing dependency: {import_name}")
        raise MissingDependency(import_name, hint)


def check_package(import_name: str) -> bool:
    """
    Check if a package is available without raising an exception.
    
    Parameters
    ----------
    import_name : str
        The Python import name of the package.
    
    Returns
    -------
    bool
        True if package can be imported, False otherwise.
    """
    try:
        importlib.import_module(import_name)
        return True
    except ImportError:
        return False


def check_dependencies(
    packages: Dict[str, Optional[str]],
    prefer_conda: bool = False,
) -> Tuple[List[str], List[str]]:
    """
    Check multiple package dependencies.
    
    Parameters
    ----------
    packages : dict
        Dictionary mapping import names to pip package names.
        Example: {'igraph': 'python-igraph', 'skmisc': 'scikit-misc'}
    prefer_conda : bool
        If True and not in venv, suggest conda for installations.
    
    Returns
    -------
    tuple of (list, list)
        (available_packages, missing_packages)
    """
    available = []
    missing = []
    
    for import_name, pip_package in packages.items():
        if check_package(import_name):
            available.append(import_name)
            logger.debug(f"Package '{import_name}' is available")
        else:
            missing.append(import_name)
            logger.warning(f"Package '{import_name}' is missing")
    
    return available, missing


def get_missing_dependencies_message(
    packages: Dict[str, Optional[str]],
    prefer_conda: bool = False,
) -> str:
    """
    Get a formatted message about missing dependencies with install instructions.
    
    Parameters
    ----------
    packages : dict
        Dictionary mapping import names to pip package names.
    prefer_conda : bool
        If True and not in venv, suggest conda installations.
    
    Returns
    -------
    str
        Formatted error message with installation instructions.
    """
    available, missing = check_dependencies(packages, prefer_conda)
    
    if not missing:
        return ""
    
    # Build error message
    lines = ["Missing required dependencies:"]
    lines.append("")
    
    for import_name in missing:
        pip_package = packages.get(import_name, import_name)
        hint = get_install_hint(import_name, pip_package, prefer_conda)
        lines.append(f"📦 {import_name}")
        for line in hint.split('\n'):
            lines.append(f"   {line}")
        lines.append("")
    
    return "\n".join(lines)
