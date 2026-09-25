# How to Contribute ✨

We warmly welcome contributions from the community! Whether you're interested in suggesting new features, fixing typos in the documentation, or making minor changes, your input is highly appreciated. 🌟

## Starting with an Issue

If you have ideas for new features or improvements, initiating a discussion in an issue is a great way to collaborate with the development team. This allows us to evaluate and discuss your suggestions together. 💡

For minor changes like fixing typos or making small edits, feel free to create a new Pull Request (PR) directly with your proposed changes. 


## Setting Up the Development Environment

1. **Fork the Repository:** 
    Start by forking the repository to your GitHub account. 🍴

2. **Clone the Forked Repository:** 
    Clone your forked repository to your local machine.

3. **Get an Environment:** Create an environment with all PPanGGOLiN prerequisites installed. For that, you can follow installation instructions [here](../user/install.md#installing-from-source-code-github).

4. **Branch from 'dev':** 
    Begin your changes from the 'dev' branch, where we incorporate changes for the upcoming release.


5. **Install in Editable Mode:** 

    To enable code editing and testing of new functionality, you can install PPanGGOLiN in editable mode using the following command:

    ```bash
    pip install -e .
    ```

    This allows you to modify the code and experiment with new features directly. 

6. **Apply Code Formatting with Black:** 
    We have integrated [Black](https://github.com/psf/black)  as our code formatter to maintain consistent styling. Code changes are automatically checked via a GitHub Action in our CI, so **ensure your code is formatted with Black before committing**.


    ```{tip}
    Integrate Black with your IDE to automatically format your changes and avoid formatting-related CI failures. 
    ```


## Making Your Changes

Keep it consistent! Match the existing code style, add docstrings to describe functions, and specify argument types.

## Update Documentation

Update docs to reflect changes—clear descriptions and examples are always helpful!

## Tests

### Running the Test Suite

We use **pytest** for both unit and functional tests.

Custom options are available to control execution:

* `--cpu=N` → number of CPUs to use in functional tests (default: 2).
* `--full` → run both unit tests and functional tests (default: unit tests only).
* `--update-golden` → update golden hashes instead of just checking them (used when expected output has legitimately changed).

Example: run the full test suite on 12 CPUs with coverage:

```bash
pytest --full -v --cpu 12
```

This will generate a coverage report in the terminal and a detailed HTML version in `htmlcov/index.html`.

### Caching of Functional Test Data

Some functional tests require building pangenomes, which is computationally expensive.
To avoid recomputing them every time, we use **pytest’s cache system**.

* On first run, the pangenome is built and stored in the pytest cache.
* On subsequent runs, the cached directory is reused, saving a lot of time.
* To force a fresh build (e.g., after code changes affecting pangenome creation), clear the cache with:

```bash
pytest --cache-clear --full --cpu 12
```

This ensures pangenomes are recomputed from scratch.


### Golden Files

Some tests compare generated files against **golden references** (checksums).
If your changes intentionally modify output formats, you can update the golden references with:

```bash
pytest --full --update-golden
```

This will rewrite the stored golden hashes with the new outputs.


### In CI

The GitHub Actions CI workflow automatically runs pytest with the appropriate options to validate new code.
If you add a new feature, make sure to also add corresponding tests under `tests/` so it gets covered locally and in CI.


## Creating a Pull Request

Once you've made your changes:

1. **Create a Pull Request:** Submit a pull request from your forked repository to the 'dev' branch on GitHub. 🚀

2. **Describe Your Changes:** Clearly describe the modifications you've made and link any associated issue(s) in the PR description. 📝

3. **Collaborative Review:** Our team will review your changes, offer feedback, and engage in discussions until we collectively agree on the implementation. 🤝

We greatly appreciate your contributions and look forward to collaborating with you! 🙌 
