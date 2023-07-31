# How to contribute to CCAKE

Thank you for your interest on contributing to CCAKE. We welcome and appreciate your interest and your effort. The guide below shows the steps one should take when contributing a feature or bug fix to the code.

## 1. Preparation

1. Whatever it is a bug fix or a feature that is to be added, please [create a new issue](https://github.com/the-nuclear-confectionery/CCAKE/issues/new).
2. From the development branch, create a new branch. This can be done locally on your local repository as
    - Go to `dev` branch: `git checkout dev`
    - Create the new branch: `git branch <branch-name>`
    - Checkout to new branch: `git checkout <branch-name>`
    - Push changes to remote: `git push origin <branch-name>`
3. Associate the branch to the issue being created

## 2. Development

1. Perform a subset of your small changes (e.g., implement a single new function)
2. Add changes to the staging area of your local repository: `git add <file/paths>`
    - Tip: You can look at the list of files changed 
3. Commit the changes: `git commit -m "TYPE: Commit description`, where `TYPE = BUG | FEAT | DOC` according to the type of changing being introduced (BUG for bugfixes, FEAT for new features and DOC for documentation being added).
4. Push changes to the feature branch with `git push origin <branch-name> `
5. Repeat this cycle until the feature development is completed

## 3. Finalize the contribution

1. Create a pull request to the `dev` branch. You can do this from within github interface
2. If the `dev` branch has progressed since the feature branch forked from it, you may need to manually solve the merge conflict.
    - In these cases, it is useful to coordinate with other developers to avoid overwriting changes.
3. Once the feature is merged on test, one should delete the feature branch in github.
    - You should also delete the branch in your local repository with `git branch -d <branch-name>`
4.  Describe what you did on the issue and close it

## 4. For maintainers only

The maintainers are the ones responsible to ensure the stability and correctness of the code. They are the code guardians. It is the maintainers responsibility to periodically look at the `dev` branch, evaluate its status and then merge into the main code.

### 4.1 Versioning

As part of the merging the `dev` into `main`, the developer should tag it with a version number. Version semantics are in the form MAJOR.MINOR.PATCH, for example `1.3.2`.

One bumps the **major** number when significant changes to the code was introduced. These typically means classes got refactored. Some may be completely gone or their functions significantly changed.

The **minor** number gets bumped when a new feature is introduced, but without large changes being made to the code.

Lastly, the **patch** number gets bumped when some small bugfixes are implemented.

As an example:
- The code is serial and a large refactoring is made to incorporate parallelism. This will result in a bump on the *major* number
- One adds a new method to compute the particlization hypersurface. This is a bump in the *minor* number.
- Lastly, one realizes that changing a given parameter in the input file is not producing any effect. This is a bug that is fixed. The *patch* number gets incremented.

Some of these tags (but not necessarily all) can be packaged to distribute a release version of the code.

For details, see [Semantic Versioning 2.0.0](https://semver.org/).