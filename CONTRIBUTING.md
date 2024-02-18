# Contributing to Our Project

We love your input! We want to make contributing to this project as easy and transparent as possible, whether it's:

- Reporting a bug
- Discussing the current state of the code
- Submitting a fix
- Proposing new features

## We Develop with GitHub

We use GitHub to host code, to track issues and feature requests, as well as accept pull requests.


## We Use [Github Flow](https://docs.github.com/en/get-started/using-github/github-flow)

So all code changes happen through pull requests. We actively welcome your pull requests!

Follow this flow to propose code changes:

1. First, check [here](https://github.com/Priusds/ParGeo/issues/), if your issue is already addressed, if so you can comment or propose code changes there.

2. If there is no related issue, [**open a new issue**](https://github.com/Priusds/ParGeo/issues/new/choose). Here you can either choose an issue template and adapt it, or open a blank one.

3. [**Create a branch**](https://docs.github.com/en/issues/tracking-your-work-with-issues/creating-a-branch-for-an-issue) from the issue.

4. Include the changes in your branch. Use the [Makefile](#the-projects-makefile) to make local code quality checks!

5. Open a **pull request** to merge your input into the master branch. Adapt the pull request to your code changes, this contains among other things an explaination of what you changed.

6. We use GitHub actions to **run code quality checks**, make sure those pass.

7. **Request a review**. The master branch is protected, everyone needs to request a review for a pull request.

8. Once the review process is done, and the reviewer approves the pull request, the proposed changes are merged into the master branch, the old branch is removed and the issue is closed.

## The Project's Makefile

Makefiles are used to automate processes. This is especially useful to test locally if your code meets the project's code requirements, meaning it passes the *ruff* linter, the *mypy* type checks and the *pytests*.

To use the Makefile you will first need an active python virual environment. You could create one by running:

```bash
python3 -m venv .venv
.venv/bin/activate
(.venv)Pargeo$ make <command> 
```

Then you can use the Makefile by running `make <command>` at the root level of the ParGeo project. Replace `<command>` with one of the following available commands:

- **pip_install:** install ParGeo using *pip*.
- **install:** install ParGeo using *poetry*.
- **format:** runs the code formatter.
- **lint:** runs the code linter and type checker.
- **test:** runs unit tests.
- **docs:** builds the documentation locally.
- **clean:** removes all files with *.geo_unrolled* and *.msh* suffix.

## Contribution License
Any contributions you make will be under the same [MIT License](LICENSE) that covers the project.
