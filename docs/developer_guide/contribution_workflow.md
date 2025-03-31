# Contribution Workflow

This document outlines the workflow for contributing to `pyremap`.

## Step 1: Fork the Repository

1. Navigate to the [pyremap repository](https://github.com/MPAS-Dev/pyremap).
2. Click the "Fork" button in the top-right corner to create your own copy of the repository.

## Step 2: Clone Your Fork

1. Clone your forked repository to your local machine:
   ```bash
   git clone https://github.com/your-username/pyremap.git
   ```
2. Navigate to the project directory:
   ```bash
   cd pyremap
   ```

## Step 3: Create a Feature Branch

1. Create a new branch for your changes:
   ```bash
   git checkout -b your-feature-name
   ```

## Step 4: Use a Worktree (Alternative to Step 3)

1. Add a new worktree for your feature branch:
   ```bash
   git worktree add ../your-feature-name
   ```
2. Navigate to the new worktree directory:
   ```bash
   cd ../your-feature-name
   ```

   > **Note:** This allows you to work on multiple branches simultaneously in separate directories.

## Step 5: Make Changes

1. Implement your changes in the codebase.
2. Write or update tests to cover your changes.
3. Update documentation if necessary.

## Step 6: Commit Your Changes

1. Stage your changes:
   ```bash
   git add .
   ```
2. Commit your changes with a descriptive message:
   ```bash
   git commit -m "Add feature: your-feature-name"
   ```

   > **Note:** Please commit often and in small pieces.  You can
   always `git rebase` later to consolidate pieces of your work into fewer
   commits, but breaking commits into smaller pieces is very difficult.

## Step 7: Push Your Changes

1. Push your feature branch to your fork:
   ```bash
   git push origin your-feature-name
   ```

## Step 8: Submit a Pull Request

1. Navigate to your fork on GitHub.
2. Click the "Compare & pull request" button.
3. Provide a detailed description of your changes.
4. Submit the pull request.

## Step 9: Address Feedback

1. Review comments from maintainers on your pull request.
2. Make necessary changes and push updates to your feature branch:
   ```bash
   git push origin your-feature-name
   ```

## Additional Resources

- Refer to the [Contribution Guidelines](index.md#contribution-guidelines) for more details.
- For help with Git, see the [Git documentation](https://git-scm.com/doc).
