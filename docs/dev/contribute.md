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

5. **Install in Editable Mode:** To enable seamless code editing and testing of new functionality, install PPanGGOLiN in editable mode using the following command:

    ```bash
    pip install -e .
    ```

    This allows you to modify the code and experiment with new features directly. 

    ```{note}
    Note: Currently, we are not utilizing any auto formatters (like autopep8 or black). Kindly refrain from using them, as it could introduce extensive changes across the project, making code review challenging for us.
    ```

## Making Your Changes

We encourage consistency in code formatting; when adding new code, try to follow the existing code structure as closely as possible. Functions should include descriptive docstrings explaining their purpose and detailing the parameters. Ensure that argument types are specified in the function definitions. 

## Update Documentation

It's essential to update the documentation to reflect your changes. Provide clear descriptions and, if necessary, examples of commands and their respective outputs.

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
