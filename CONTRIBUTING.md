# Contributing to Novae-Seurat-GUI

Thank you for your interest in contributing to Novae-Seurat-GUI! This document provides guidelines for contributing to the project.

## Code of Conduct

By participating in this project, you agree to maintain a respectful and inclusive environment for all contributors.

## How to Contribute

### Reporting Bugs

If you find a bug, please open an issue with:
- A clear description of the problem
- Steps to reproduce
- Expected vs actual behavior
- Your environment (Python version, OS, etc.)

### Suggesting Enhancements

We welcome suggestions for new features or improvements. Please open an issue with:
- A clear description of the enhancement
- Use cases and benefits
- Any implementation ideas (optional)

### Pull Requests

1. **Fork the repository** and create your branch from `main`

2. **Set up development environment**:
   ```bash
   git clone https://github.com/yourusername/XSpatialNovae.git
   cd XSpatialNovae
   pip install -e ".[dev]"
   ```

3. **Make your changes**:
   - Write clear, documented code
   - Follow existing code style
   - Add tests for new functionality
   - Update documentation as needed

4. **Test your changes**:
   ```bash
   # Run tests
   pytest tests/
   
   # Check code formatting
   black novae_seurat_gui tests
   
   # Check linting
   flake8 novae_seurat_gui
   ```

5. **Commit your changes**:
   - Use clear, descriptive commit messages
   - Reference issue numbers if applicable
   ```bash
   git commit -m "Add feature X to address #123"
   ```

6. **Push to your fork** and submit a pull request

7. **Wait for review**:
   - Maintainers will review your PR
   - Address any requested changes
   - Once approved, your PR will be merged

## Development Guidelines

### Code Style

- Follow PEP 8 guidelines
- Use Black for code formatting (line length: 100)
- Use type hints where appropriate
- Write docstrings for all public functions/classes

### Testing

- Write unit tests for new functions
- Maintain or improve code coverage
- Test edge cases and error handling
- Use pytest fixtures for common setup

### Documentation

- Update README.md for user-facing changes
- Add docstrings with clear parameter descriptions
- Include examples in docstrings when helpful
- Update configuration file documentation

### Commit Messages

Use conventional commit format:
- `feat:` for new features
- `fix:` for bug fixes
- `docs:` for documentation changes
- `test:` for adding/updating tests
- `refactor:` for code refactoring
- `style:` for formatting changes
- `chore:` for maintenance tasks

Example: `feat: add proteomics preprocessing pipeline`

## Project Structure

```
XSpatialNovae/
├── novae_seurat_gui/     # Main package
│   ├── io/               # Data loading and validation
│   ├── qc/               # Quality control
│   ├── spatial/          # Spatial analysis
│   ├── modeling/         # Novae model wrappers
│   ├── niche/            # Niche analysis
│   ├── viz/              # Visualization
│   ├── export/           # R-friendly exports
│   └── cli.py            # Command-line interface
├── app.py                # Streamlit GUI
├── configs/              # Configuration templates
├── examples/             # Example scripts
├── tests/                # Unit tests
└── docs/                 # Documentation (future)
```

## Running Tests Locally

```bash
# Run all tests
pytest

# Run specific test file
pytest tests/test_io.py

# Run with coverage
pytest --cov=novae_seurat_gui --cov-report=html

# View coverage report
open htmlcov/index.html
```

## Building Documentation

(To be added when documentation is set up)

## Questions?

- Open an issue for questions
- Check existing issues and PRs
- Reach out to maintainers if needed

Thank you for contributing!
