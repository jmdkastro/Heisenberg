# Contributing to Heisenberg

Thank you for your interest in contributing to Heisenberg! This document provides guidelines and instructions for contributors.

## Getting Started

### Development Environment Setup

1. **Fork and clone the repository**

```bash
git clone https://github.com/mustang-project/heisenberg.git
cd heisenberg
```

2. **Create a virtual environment**

```bash
python -m venv venv
source venv/bin/activate  # Linux/macOS
# or: venv\Scripts\activate  # Windows
```

3. **Install development dependencies**

```bash
cd Python/
pip install -e ".[dev]"
```

4. **Verify the installation**

```bash
pytest tests/
```

## Code Style

### Python

- **Style Guide**: Follow [PEP 8](https://www.python.org/dev/peps/pep-0008/)
- **Formatting**: Use [black](https://github.com/psf/black) with default settings
- **Import Sorting**: Use [isort](https://pycqa.github.io/isort/) with black-compatible profile
- **Type Hints**: Use type hints for all public functions
- **Docstrings**: Use [NumPy style](https://numpydoc.readthedocs.io/en/latest/format.html) docstrings

Format your code before committing:

```bash
black heisenberg/
isort heisenberg/
```

### IDL

- Follow existing code style in the repository
- Use 3-space indentation (matching existing files)
- Include header comments for all procedures and functions
- Document parameters in comments

### Docstring Example

```python
def fit_tuningfork(
    flux_ratios: np.ndarray,
    apertures: np.ndarray,
    tstar: float,
    config: Optional[FitConfig] = None,
) -> FitResult:
    """
    Fit the KL14 uncertainty principle model to flux ratio data.

    Parameters
    ----------
    flux_ratios : np.ndarray
        Observed flux ratios (stellar/gas) for each aperture size.
        Shape: (n_apertures,)
    apertures : np.ndarray
        Aperture diameters in parsecs. Shape: (n_apertures,)
    tstar : float
        Reference stellar tracer timescale in Myr.
    config : FitConfig, optional
        Fitting configuration. If None, uses defaults.

    Returns
    -------
    FitResult
        Fitted parameters including tgas, tover, lambda with errors.

    Raises
    ------
    ValueError
        If flux_ratios and apertures have different lengths.
    FitConvergenceError
        If the fitting algorithm fails to converge.

    Examples
    --------
    >>> ratios = np.array([0.5, 0.6, 0.7, 0.8, 0.9])
    >>> apertures = np.array([50, 100, 200, 400, 800])
    >>> result = fit_tuningfork(ratios, apertures, tstar=4.0)
    >>> print(f"t_gas = {result.tgas:.1f} Myr")
    t_gas = 10.5 Myr
    """
```

## Testing

### Running Tests

```bash
# Run all tests
pytest tests/

# Run with coverage
pytest tests/ --cov=heisenberg --cov-report=html

# Run specific test file
pytest tests/unit/test_fitting.py

# Run specific test
pytest tests/unit/test_fitting.py::test_fit_converges

# Run with verbose output
pytest tests/ -v
```

### Writing Tests

- Place unit tests in `tests/unit/`
- Place integration tests in `tests/integration/`
- Use descriptive test names: `test_<function>_<scenario>_<expected_result>`
- Use pytest fixtures for common setup
- Aim for >80% coverage on new code

Test example:

```python
import pytest
import numpy as np
from heisenberg.core.fitting import fit_tuningfork, FitResult

class TestFitTuningfork:
    """Tests for the fit_tuningfork function."""

    def test_returns_fit_result(self):
        """fit_tuningfork should return a FitResult object."""
        ratios = np.array([0.5, 0.6, 0.7, 0.8, 0.9])
        apertures = np.array([50, 100, 200, 400, 800])

        result = fit_tuningfork(ratios, apertures, tstar=4.0)

        assert isinstance(result, FitResult)

    def test_tgas_positive(self):
        """Fitted tgas should be positive."""
        ratios = np.array([0.5, 0.6, 0.7, 0.8, 0.9])
        apertures = np.array([50, 100, 200, 400, 800])

        result = fit_tuningfork(ratios, apertures, tstar=4.0)

        assert result.tgas > 0

    def test_mismatched_lengths_raises(self):
        """Should raise ValueError for mismatched array lengths."""
        ratios = np.array([0.5, 0.6, 0.7])
        apertures = np.array([50, 100])  # Different length

        with pytest.raises(ValueError):
            fit_tuningfork(ratios, apertures, tstar=4.0)
```

## Pull Request Process

1. **Create a feature branch**

```bash
git checkout -b feature/my-new-feature
```

2. **Make your changes**
   - Write code following the style guidelines
   - Add or update tests
   - Update documentation if needed

3. **Run tests and checks**

```bash
# Run tests
pytest tests/

# Check formatting
black --check heisenberg/
isort --check-only heisenberg/

# Run linter
flake8 heisenberg/
```

4. **Commit your changes**

```bash
git add .
git commit -m "Add feature: description of changes"
```

Use clear, descriptive commit messages:
- Start with a verb (Add, Fix, Update, Remove, Refactor)
- Keep the first line under 72 characters
- Add details in the body if needed

5. **Push and create PR**

```bash
git push origin feature/my-new-feature
```

Then create a Pull Request on GitHub with:
- Clear title describing the change
- Description of what and why
- Reference any related issues

6. **Address review feedback**

Make requested changes and push updates to the same branch.

## Types of Contributions

### Bug Reports

- Use the GitHub issue tracker
- Include: Python/IDL version, OS, error message, steps to reproduce
- Minimal reproducible example if possible

### Feature Requests

- Open a GitHub issue describing the feature
- Explain the use case and expected behavior
- Discuss before implementing large features

### Documentation

- Fix typos or unclear explanations
- Add examples
- Improve API documentation
- Translate documentation

### Code Contributions

- Bug fixes
- New features (discuss first for large changes)
- Performance improvements
- Test coverage improvements

## License

By contributing to Heisenberg, you agree that your contributions will be licensed under the **GNU General Public License v3.0** (GPL-3.0), the same license as the rest of the project.

See [LICENSE](LICENSE) for the full license text.

## Questions?

- Open a GitHub Discussion for questions
- Check existing issues and discussions first
- Be respectful and constructive

Thank you for contributing to Heisenberg!
