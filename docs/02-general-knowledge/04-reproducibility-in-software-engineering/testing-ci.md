---
title: Testing & Continuous Integration
slug: /general-knowledge/testing-ci
sidebar_position: 5
---

:::caution Starting point, contributions welcome
This page is deliberately short and covers only the essentials. If you have experience with testing research software and would like to extend it, please [contact us](mailto:rdhub@mdsi.tum.de) or use the "Edit this page" link below.
:::

# **Testing & Continuous Integration**

Research results are only as reliable as the code that produced them. A small error in an analysis script does not crash; it silently produces wrong numbers. **Tests** are small programs that check that your code does what you expect. **Continuous integration (CI)** runs them automatically every time the code changes.

You do not need full test coverage to benefit. A handful of well-chosen tests catches most accidents, and makes it safe to change code later.

## **1\. What to Test First**

| Kind of test | What it checks | A good first step |
|---|---|---|
| **Unit test** | One function gives the right output for a known input | Test the functions that implement your method or do calculations, including edge cases (empty input, missing values, zeros) |
| **End-to-end test** | The whole analysis or workflow runs on a small example dataset | Keep a tiny test dataset in the repository that runs in seconds to minutes |
| **Regression test** | Results have not changed unintentionally | Store the expected output for the test dataset and compare against it |

Practical rules:

- **Start with the code that matters most** for your conclusions, not with the code that is easiest to test.
- **Turn every bug into a test.** When you fix an error, first write a test that fails because of it.
- **Compare numbers with a tolerance,** not for exact equality, and set random seeds so that results are repeatable.
- **Keep tests fast,** otherwise nobody runs them.

## **2\. Tools**

Use the standard framework of your language: [pytest](https://docs.pytest.org/) for Python, [testthat](https://testthat.r-lib.org/) for R, the built-in `Test` module for Julia, GoogleTest or Catch2 for C++. A minimal Python example:

```python
# test_stats.py
import pytest
from mypackage.stats import normalise

def test_normalise_sums_to_one():
    result = normalise([1.0, 1.0, 2.0])
    assert sum(result) == pytest.approx(1.0)

def test_normalise_rejects_empty_input():
    with pytest.raises(ValueError):
        normalise([])
```

Running `pytest` in the project folder finds and executes all such tests.

For [workflows](/general-knowledge/wms), test at the workflow level as well: a dry run (`snakemake -n`) checks that the workflow is consistent, and a run on a small test dataset checks that it works. Community pipelines such as nf-core ship a `test` profile for exactly this purpose.

## **3\. Continuous Integration**

A CI service runs a defined set of commands on a fresh machine whenever you push changes or open a merge request. This also answers the question "does it work on a computer other than mine?".

- **GitHub** provides [GitHub Actions](https://docs.github.com/en/actions), free for public repositories.
- **GitLab**, including [LRZ GitLab](https://gitlab.lrz.de), provides [GitLab CI/CD](https://docs.gitlab.com/ci/). Jobs need a *runner* to execute on; check the LRZ documentation for which runners are available to your project.

A minimal GitHub Actions configuration for a Python project, saved as `.github/workflows/tests.yml`:

```yaml
name: tests
on: [push, pull_request]

jobs:
  test:
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v4
      - uses: astral-sh/setup-uv@v6
      - run: uv sync --locked
      - run: uv run pytest
```

Installing from the lock file (see [Environments and Dependency Management](/general-knowledge/containers#first-step-environments-and-dependency-management)) makes sure the tests run with the same package versions you use.

CI is not limited to tests. Typical further jobs are code formatting and static checks ("linting"), building [container](/general-knowledge/containers) images, and building documentation. This Knowledge Hub is an example: every change is built automatically, and the build fails if a page contains a broken internal link.

:::tip Do not run tests on sensitive data in public CI
CI jobs on public services run on external servers. Use synthetic or openly licensed test data, and never store passwords or access tokens in the repository; use the secret store of the CI service instead.
:::

## **Further Information**

- [The Turing Way: Code Testing](https://book.the-turing-way.org/reproducible-research/testing/) and [Continuous Integration](https://book.the-turing-way.org/reproducible-research/ci/)
- Wilson G. et al. (2017): [Good enough practices in scientific computing](https://doi.org/10.1371/journal.pcbi.1005510). *PLoS Computational Biology*.
- [CodeRefinery: Software testing](https://coderefinery.github.io/testing/) (hands-on lesson)
