# Releasing a New Version

This document describes the steps for maintainers to tag and release a new
version of `pyremap`, and to update the conda-forge feedstock.

## Version Bump and Dependency Updates

1. **Update the Version Number**

   - Manually update the version number in the following files:

     - `pyremap/version.py`
     - `ci/recipe/meta.yaml`

   - Make sure the version follows [semantic versioning](https://semver.org/).
     For release candidates, use versions like `2.1.0rc1`.

2. **Check and Update Dependencies**

   - Ensure that dependencies and their constraints are up-to-date and
     consistent in:

     - `ci/recipe/meta.yaml` (for the conda-forge release)
     - `pyproject.toml` (for PyPI; used for sanity checks)
     - `dev-spec.txt` (development dependencies; should be a superset)

   - Use the GitHub "Compare" feature to check for dependency changes between
     releases:
     [https://github.com/MPAS-Dev/pyremap/compare](https://github.com/MPAS-Dev/pyremap/compare)

3. **Make a PR and merge it**

   - Open a PR for the version bump and dependency changes and merge once
     approved.

## Tagging and Publishing a Release Candidate

4. **Tagging a Release Candidate**

   - For release candidates, **do not create a GitHub release page**. Just
     create a tag from the command line:

     - Make sure your changes are merged into `main` or your own update branch
       (e.g. `update-to-2.1.0`) and your local repo is up to date.

     - Tag the release candidate (e.g., `2.1.0rc1`):

       ```bash
       git checkout main
       git fetch --all -p
       git reset --hard origin/main
       git tag 2.1.0rc1
       git push origin 2.1.0rc1
       ```

       (Replace `2.1.0rc1` with your actual version, and `main` with your
       branch if needed.)

     **Note:** This will only create a tag. No release page will be created on
     GitHub.

## Updating the conda-forge Feedstock for a Release Candidate

5. **Manual Feedstock Update (Required for Release Candidates)**

   - The conda-forge feedstock does **not** update automatically for release
     candidates.
   - You must always create a PR manually, and it must target the `dev`
     branch of the feedstock.

   Steps:

   - Download the release tarball:

     ```bash
     wget https://github.com/MPAS-Dev/pyremap/archive/refs/tags/<version>.tar.gz
     ```

   - Compute the SHA256 checksum:

     ```bash
     shasum -a 256 <version>.tar.gz
     ```

   - In the `meta.yaml` of the feedstock recipe:
     - Set `{% set version = "<version>" %}`
     - Set the new `sha256` value
     - Update dependencies if needed

   - Commit, push to a new branch, and open a PR **against the `dev` branch**
     of the feedstock:
     https://github.com/conda-forge/pyremap-feedstock

   - Follow any instructions in the PR template and merge once approved

## Releasing a Stable Version

6. **Publishing a Stable Release**

   - For stable releases, create a GitHub release page as follows:

     - Go to [https://github.com/MPAS-Dev/pyremap/releases](https://github.com/MPAS-Dev/pyremap/releases)
     - Click "Draft a new release"
     - Enter a tag (e.g., `2.1.0`)
     - Set the release title to the version prefixed with `v` (e.g., `v2.1.0`)
     - Generate or manually write release notes
     - Click "Publish release"

## Updating the conda-forge Feedstock for a Stable Release

7. **Automatic Feedstock Update (Preferred Method)**

   - Wait for the `regro-cf-autotick-bot` to open a PR at:
     [https://github.com/conda-forge/pyremap-feedstock](https://github.com/conda-forge/pyremap-feedstock)

   - This may take several hours to a day.

   - Review the PR:
     - Confirm the version bump and dependency changes
     - Merge once CI checks pass

   **Note:** If you are impatient, you can accelerate this process by creating
   a bot issue at:
   [https://github.com/conda-forge/pyremap-feedstock/issues](https://github.com/conda-forge/pyremap-feedstock/issues)
   with the subject `@conda-forge-admin, please update version`. This will open
   a new PR with the version within a few minutes.

   - If the bot PR does not appear or is too slow, you may update manually (see
     the manual steps for release candidates above, but target the `main`
     branch of the feedstock).

## Post Release Actions

8. **Verify and Announce**

   - Install the package in a clean environment to test:

     ```bash
     conda create -n test-pyremap -c conda-forge pyremap=<version>
     ```

   - Optionally announce the release on relevant communication channels

   - Update any documentation or release notes as needed
