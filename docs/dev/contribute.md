# How to Contribute ✨

We warmly welcome contributions from the community! Whether you're interested in suggesting new features, fixing typos in the documentation, or making minor changes, your input is highly appreciated. 🌟

## Starting with an Issue

If you have ideas for new features or improvements, initiating a discussion in an issue is a great way to collaborate with the development team. This allows us to evaluate and discuss your suggestions together. 💡

For minor changes like fixing typos or making small edits, feel free to create a new Pull Request (PR) directly with your proposed changes. 


## Setting Up the Development Environment

1. **Fork the Repository:** Start by forking the repository to your GitHub account. 🍴

2. **Clone the Forked Repository:** Clone your forked repository to your local machine.

3. **Get an Environment:** Create an environment with all PPanGGOLiN prerequisites installed. For that, you can follow installation instructions [here](../user/install.md#installing-from-source-code-github).

4. **Branch from 'dev':** Begin your changes from the 'dev' branch, where we incorporate changes for the upcoming release.

5. **Install in Editable Mode:** To enable code editing and testing of new functionality, you can install PPanGGOLiN in editable mode using the following command:

    ```bash
    pip install -e .
    ```

    This allows you to modify the code and experiment with new features directly. 

6. **Apply Code Formatting with Black:** We have integrated Black as our code formatter to maintain consistent styling. Code changes are automatically checked via a GitHub Action in our CI, so ensure your code is formatted with Black before committing.


    ```{tip}
    Integrate Black with your IDE to automatically format your changes and avoid formatting-related CI failures. 
    ```


## Making Your Changes

Keep it consistent! Match the existing code style, add docstrings to describe functions, and specify argument types.

## Update Documentation

Update docs to reflect changes—clear descriptions and examples are always helpful!

## Tests

### Continuous Integration (CI) Workflow

We've set up a CI workflow in the Actions tab that executes a series of PPanGGOLiN commands to validate their functionality, and also compares the contents of the PPanGGOLiN info files generated during the workflow with the expected ones stored in the `testingDataset` directory. If you've introduced a new feature, consider adding a command line to the CI YAML file to test it and ensure its seamless integration.

The CI workflow can be launched locally using the Python script `launch_test_locally.py` located in the `testingDataset` directory. This script reads the CI pipeline file and creates a bash script to facilitate local execution of the pipeline.


To setup the local execution with 10 CPUs in the local_CI directory, execute the following command:

```bash
python testingDataset/launch_test_locally.py -o local_CI -c 10 

```

Then, run the local CI using the following command:

```bash
(cd local_CI; bash launch_test_command.sh)
```

### Unit Tests

While not mandatory for all PPanGGOLiN code, incorporating unit tests for your additions can be advantageous. The test suite is located in the 'tests' directory at the root of the project.

## Creating a Pull Request

Once you've made your changes:

1. **Create a Pull Request:** Submit a pull request from your forked repository to the 'dev' branch on GitHub. 🚀

2. **Describe Your Changes:** Clearly describe the modifications you've made and link any associated issue(s) in the PR description. 📝

3. **Collaborative Review:** Our team will review your changes, offer feedback, and engage in discussions until we collectively agree on the implementation. 🤝

We greatly appreciate your contributions and look forward to collaborating with you! 🙌 
